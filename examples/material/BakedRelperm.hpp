// Baked-table relperm/pc evaluation (experiment; see relperm_rewrite_plan.md).
//
// For each cell the true material law (multiplexer params from the
// EclMaterialLawManager, i.e. EPS included) is sampled onto minimal exact
// piecewise-linear tables on three axes:
//   Sw axis: [krw, pcow-diff]     (pcow-diff = pC[w]-pC[o])
//   Sg axis: [krg, pcgo-diff]
//   S  axis: [kro_ow, kro_go]     (S = max(Sw,Swco)+Sg; ECL default 3ph oil)
// Exactness comes from bake-then-refine: candidate nodes are densified until
// midpoint probes of every interval match the true law, then collinear nodes
// are pruned. Tables are deduplicated across cells by (satnum, scaled EPS
// info), so storage is a pool + one index per cell — a GPU-friendly SoA shape.
#ifndef OPM_BENCH_BAKED_RELPERM_HPP
#define OPM_BENCH_BAKED_RELPERM_HPP

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidmatrixinteractions/EclStone1Material.hpp>
#include <opm/material/fluidmatrixinteractions/EclStone2Material.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm::Bench {

// Minimal fluid state carrying only saturations (all the family-I curves read).
template <class Scalar, int numPhases>
struct SatOnlyFluidState
{
    std::array<Scalar, numPhases> sat{};
    Scalar saturation(unsigned phaseIdx) const
    { return sat[phaseIdx]; }
};

// One two-column curve with precomputed slopes; x strictly increasing on [0,1].
// Node values interleaved as {v0, s0, v1, s1} so an evaluation after the
// search touches a single cache line (the flowdiagnostics table layout).
struct BakedCurve
{
    std::vector<double> x;
    // 2*nCols doubles per node: v0,s0, v1,s1[, v2,s2]. Interleaved so one
    // cache line after the search yields every column and its slope.
    std::vector<double> vals;
    int nCols = 2;

    // exact search acceleration: a uniform hint grid maps the abscissa to a
    // starting node; the forward scan is bounded by nodes per hint cell.
    static constexpr std::size_t hintCells = 64;
    double x0 = 0.0;
    double invH = 0.0;
    std::array<std::uint16_t, hintCells + 1> hint{};

    void buildHint()
    {
        x0 = x.front();
        const double span = x.back() - x0;
        invH = (span > 0.0) ? double(hintCells) / span : 0.0;
        std::size_t i = 0;
        for (std::size_t j = 0; j <= hintCells; ++j) {
            const double xj = x0 + (invH > 0.0 ? double(j) / invH : 0.0);
            while (i + 2 < x.size() && x[i + 1] <= xj) {
                ++i;
            }
            hint[j] = static_cast<std::uint16_t>(i);
        }
    }

    // one hint lookup + short scan serves every column; derivatives analytic
    template <class Evaluation>
    void evalCols(const Evaluation& xe, Evaluation* out, int n) const
    {
        double xd;
        if constexpr (std::is_floating_point_v<Evaluation>) {
            xd = xe;
        }
        else {
            xd = xe.value();
        }
        const double* xb = x.data();
        const double cell = std::min(double(hintCells), std::max(0.0, (xd - x0) * invH));
        std::size_t lo = hint[static_cast<std::size_t>(cell)];
        while (lo + 2 < x.size() && xb[lo + 1] <= xd) {
            ++lo;
        }
        const double dx = xd - xb[lo];
        const double* nd = vals.data() + 2 * nCols * lo;
        for (int c = 0; c < n; ++c) {
            const double v = nd[2 * c];
            const double sl = nd[2 * c + 1];
            if constexpr (std::is_floating_point_v<Evaluation>) {
                out[c] = v + sl * dx;
            }
            else {
                out[c].setValue(v + sl * dx);
                for (int k = 0; k < xe.size(); ++k) {
                    out[c].setDerivative(k, sl * xe.derivative(k));
                }
            }
        }
    }

    template <class Evaluation>
    void eval2(const Evaluation& xe, Evaluation& c0, Evaluation& c1) const
    {
        double xd;
        if constexpr (std::is_floating_point_v<Evaluation>) {
            xd = xe;
        }
        else {
            xd = xe.value();
        }
        const double* xb = x.data();
        const double cell = std::min(double(hintCells), std::max(0.0, (xd - x0) * invH));
        std::size_t lo = hint[static_cast<std::size_t>(cell)];
        while (lo + 2 < x.size() && xb[lo + 1] <= xd) {
            ++lo;
        }
        const double dx = xd - xb[lo];
        const double* nd = vals.data() + 4 * lo;
        if constexpr (std::is_floating_point_v<Evaluation>) {
            c0 = nd[0] + nd[1] * dx;
            c1 = nd[2] + nd[3] * dx;
        }
        else {
            c0.setValue(nd[0] + nd[1] * dx);
            c1.setValue(nd[2] + nd[3] * dx);
            for (int k = 0; k < xe.size(); ++k) {
                c0.setDerivative(k, nd[1] * xe.derivative(k));
                c1.setDerivative(k, nd[3] * xe.derivative(k));
            }
        }
    }
};

// Experiment (see baked_satfunc_implementation.md §7.2): the same curve
// re-tabulated on a *uniform* abscissa, so evaluation needs no search at all
// — one multiply and a truncation. Node count is chosen so the maximum
// deviation from the exact curve stays under a tolerance; that makes it a
// physically equivalent redefinition rather than an exact representation.
struct UniformCurve
{
    double x0 = 0.0;
    double invH = 0.0;
    std::size_t n = 0;
    std::vector<double> vals;   // 4 doubles per node: v0, s0, v1, s1

    template <class Evaluation>
    void eval2(const Evaluation& xe, Evaluation& c0, Evaluation& c1) const
    {
        double xd;
        if constexpr (std::is_floating_point_v<Evaluation>) {
            xd = xe;
        }
        else {
            xd = xe.value();
        }
        const double t = (xd - x0) * invH;
        const double tc = std::min(double(n - 2), std::max(0.0, t));
        const std::size_t i = static_cast<std::size_t>(tc);
        const double dx = xd - (x0 + double(i) / invH);
        const double* nd = vals.data() + 4 * i;
        if constexpr (std::is_floating_point_v<Evaluation>) {
            c0 = nd[0] + nd[1] * dx;
            c1 = nd[2] + nd[3] * dx;
        }
        else {
            c0.setValue(nd[0] + nd[1] * dx);
            c1.setValue(nd[2] + nd[3] * dx);
            for (int k = 0; k < xe.size(); ++k) {
                c0.setDerivative(k, nd[1] * xe.derivative(k));
                c1.setDerivative(k, nd[3] * xe.derivative(k));
            }
        }
    }
};

struct BakedTable
{
    BakedCurve sw;   // krw, pcow-diff
    BakedCurve sg;   // krg, pcgo-diff
    BakedCurve so;   // kro_ow, kro_go   (abscissa is S = sw+sg)
    // Carlson hysteresis: only the two non-wetting relperms are hysteretic,
    // so a second pool holds just their imbibition variants (column 0 used).
    BakedCurve sgImb; // krg  from the gas-oil imbibition curve, vs sg
    BakedCurve soImb; // kro_ow from the oil-water imbibition curve, vs S
    bool hasImb = false;
    UniformCurve uSw, uSg, uSo;
    double swco = 0.0;
    // 0 = ECL default (three axes); 1 = Stone 1, 2 = Stone 2 (two axes,
    // three columns each, oil constituents in column 2)
    int approach = 0;
};

template <class MaterialLaw, class Manager, class Scalar>
class BakedRelperm
{
    using MultiplexerParams = typename Manager::MaterialLawParams;
    using DefaultMaterial = typename MaterialLaw::DefaultMaterial;
    using DefaultParams = typename DefaultMaterial::Params;
    static constexpr int numPhases = 3;
    static constexpr unsigned waterIdx = MaterialLaw::Traits::wettingPhaseIdx;
    static constexpr unsigned oilIdx = MaterialLaw::Traits::nonWettingPhaseIdx;
    static constexpr unsigned gasIdx = MaterialLaw::Traits::gasPhaseIdx;

    using FS = SatOnlyFluidState<Scalar, numPhases>;
    static constexpr bool hasHysteresis = MaterialLaw::Traits::enableHysteresis;

public:
    // truth functions of one variable each (see header comment)
    struct Truth
    {
        const DefaultParams* p;
        double swco;

        std::array<double, 2> swAxis(double swv) const
        {   // krw, pcow-diff at (Sw=swv, Sg=0)
            FS fs; fs.sat[waterIdx] = swv; fs.sat[oilIdx] = 1.0 - swv; fs.sat[gasIdx] = 0.0;
            std::array<double, numPhases> pC{};
            DefaultMaterial::capillaryPressures(pC, *p, fs);
            const double krw = DefaultMaterial::template krw<FS, double>(*p, fs);
            return { krw, pC[waterIdx] - pC[oilIdx] };
        }
        std::array<double, 2> sgAxis(double sgv) const
        {   // krg, pcgo-diff at (Sw=swco, Sg=sgv)
            FS fs; fs.sat[waterIdx] = swco; fs.sat[gasIdx] = sgv;
            fs.sat[oilIdx] = 1.0 - swco - sgv;
            std::array<double, numPhases> pC{};
            DefaultMaterial::capillaryPressures(pC, *p, fs);
            const double krg = DefaultMaterial::template krg<FS, double>(*p, fs);
            return { krg, pC[gasIdx] - pC[oilIdx] };
        }
        std::array<double, 2> soAxis(double S) const
        {   // kro_ow(S) at (Sw=S, Sg=0); kro_go(S) at (Sw=swco, Sg=S-swco)
            FS a; a.sat[waterIdx] = S; a.sat[oilIdx] = 1.0 - S; a.sat[gasIdx] = 0.0;
            FS b; b.sat[waterIdx] = swco; b.sat[gasIdx] = std::max(S - swco, 0.0);
            b.sat[oilIdx] = 1.0 - swco - b.sat[gasIdx];
            const double kroOw = DefaultMaterial::template relpermOilInOilWaterSystem<double, FS>(*p, a);
            const double kroGo = DefaultMaterial::template relpermOilInOilGasSystem<double, FS>(*p, b);
            return { kroOw, kroGo };
        }

        // Stone 1 / Stone 2 put both oil constituents on the sw and sg axes
        // (kro_ow at sw, kro_go at 1-Swco-sg), so those axes carry three
        // columns and no S axis is needed.
        // Imbibition variants of the two hysteretic curves, taken straight
        // from the imbibition parameter set (the composed law would return
        // the drainage branch while the cell state is at its initial value).
        // Arguments are the two-phase wetting saturations the ECL default
        // material passes: 1-Swl-sg for gas-oil, S for oil-water.
        std::array<double, 2> sgAxisImb(double sgv) const
        {
            if constexpr (hasHysteresis) {
                using GasOilLaw = typename DefaultMaterial::GasOilMaterialLaw;
                using EffLaw = typename GasOilLaw::EffectiveLaw;
                const double arg = 1.0 - swco - sgv;
                return { EffLaw::template twoPhaseSatKrn<double>(p->gasOilParams().imbibitionParams(), arg),
                         0.0 };
            }
            else {
                return { 0.0, 0.0 };
            }
        }
        std::array<double, 2> soAxisImb(double S) const
        {
            if constexpr (hasHysteresis) {
                using OilWaterLaw = typename DefaultMaterial::OilWaterMaterialLaw;
                using EffLaw = typename OilWaterLaw::EffectiveLaw;
                return { EffLaw::template twoPhaseSatKrn<double>(p->oilWaterParams().imbibitionParams(), S),
                         0.0 };
            }
            else {
                return { 0.0, 0.0 };
            }
        }
    };

    // per-cell factor slots: [krw, kro_ow, pcow, krg, kro_go, pcgo]
    static constexpr int numFactors = 6;

    // Stone 1 / Stone 2 put both oil constituents on the sw and sg axes
    // (kro_ow at sw, kro_go at 1-Swco-sg), so those axes carry three columns
    // and no S axis is needed.
    template <class StoneMat, class StoneParams>
    static std::array<double, 3> swAxisStone(const StoneParams& sp, double swv)
    {
        FS fs; fs.sat[waterIdx] = swv; fs.sat[oilIdx] = 1.0 - swv; fs.sat[gasIdx] = 0.0;
        std::array<double, numPhases> pC{};
        StoneMat::capillaryPressures(pC, sp, fs);
        const double krw = StoneMat::template krw<FS, double>(sp, fs);
        const double kroOw = StoneMat::template relpermOilInOilWaterSystem<double, FS>(sp, fs);
        return { krw, pC[waterIdx] - pC[oilIdx], kroOw };
    }
    template <class StoneMat, class StoneParams>
    static std::array<double, 3> sgAxisStone(const StoneParams& sp, double swco, double sgv)
    {
        FS fs; fs.sat[waterIdx] = swco; fs.sat[gasIdx] = sgv;
        fs.sat[oilIdx] = 1.0 - swco - sgv;
        std::array<double, numPhases> pC{};
        StoneMat::capillaryPressures(pC, sp, fs);
        const double krg = StoneMat::template krg<FS, double>(sp, fs);
        const double kroGo = StoneMat::template relpermOilInOilGasSystem<double, FS>(sp, fs);
        return { krg, pC[gasIdx] - pC[oilIdx], kroGo };
    }

    // Empty when the deck's saturation-function configuration is inside what
    // the baked representation reproduces exactly; otherwise a reason. A
    // caller must consult this and fall back to the generic law — the baked
    // tables silently ignore anything listed here.
    static std::string unsupportedReason(const Manager& manager)
    {
        switch (manager.threePhaseApproach()) {
        case Opm::EclMultiplexerApproach::Default: break;
        case Opm::EclMultiplexerApproach::Stone1:  return "Stone 1 three-phase kro";
        case Opm::EclMultiplexerApproach::Stone2:  return "Stone 2 three-phase kro";
        case Opm::EclMultiplexerApproach::TwoPhase: return "two-phase system";
        case Opm::EclMultiplexerApproach::OnePhase: return "single-phase system";
        }
        if (manager.hasDirectionalRelperms()) {
            return "directional relperms (KRnX/KRnY/KRnZ)";
        }
        if (manager.hasDirectionalImbnum()) {
            return "directional IMBNUM";
        }
        if (manager.enableHysteresis()) {
            if constexpr (!hasHysteresis) {
                // Traits compiled without hysteresis: the deck asks for it and
                // the baked tables would silently drop it. Worse than slow.
                return "deck enables EHYSTR but the traits disable hysteresis";
            }
            const auto& hc = manager.hysteresisConfig();
            if (hc.krHysteresisModel() != 0) {
                // 0 = Carlson with drainage krw, the only model baked. 1 is
                // Carlson but takes krw from the imbibition curve; 2-4 are
                // Killough, which rescales rather than shifts.
                return "krHysteresisModel " + std::to_string(hc.krHysteresisModel())
                     + " (only 0, Carlson, is baked)";
            }
            if (hc.pcHysteresisModel() >= 0) {
                return "capillary-pressure hysteresis";
            }
            if (hc.enableWagHysteresis()) {
                return "WAG hysteresis";
            }
        }
        return {};
    }

    // force=true builds anyway, for benchmarks that deliberately price a
    // configuration the tables do not reproduce. Never set it in a simulator.
    void build(const Manager& manager, std::size_t nc, bool force = false)
    {
        const std::string why = unsupportedReason(manager);
        if (!why.empty() && !force) {
            throw std::runtime_error("BakedRelperm cannot represent this deck: " + why);
        }
        cellTable_.resize(nc);
        cellSwco_.resize(nc);
        cellFactors_.assign(numFactors * nc, 1.0);
        if constexpr (hasHysteresis) {
            cellHyst_.assign(4 * nc, 0.0);
        }
        std::map<std::vector<std::int64_t>, std::uint32_t> dedup;
        std::vector<std::array<double, numFactors>> tableRefs;

        for (std::size_t c = 0; c < nc; ++c) {
            const auto& info = manager.oilWaterScaledEpsInfoDrainage(c);
            const int satnum = manager.satnumRegionIdx(c);
            // Vertical endpoint scales (maxKr*, maxPc*) are NOT part of the
            // key: 2pt vertical scaling — and 3pt with equal Kr*r/maxKr*
            // ratios, which the key preserves — is linear in the ordinate and
            // factored into per-cell multipliers. Without this, per-cell pc
            // (SWATINIT) or kr scaling degenerates to one table per cell.
            auto key = makeKey(satnum, info,
                               manager.oilWaterConfig(), manager.gasOilConfig());
            if constexpr (hasHysteresis) {
                // The imbibition curves carry their own endpoint scaling
                // (ISGCR, ISWCR, ... — Norne sets ISGCR), which the drainage
                // info above does not see. Without this two cells with equal
                // drainage scaling but different imbibition scaling would
                // share a table and one of them would get the other's
                // imbibition curve.
                const DefaultParams& dp =
                    manager.materialLawParams(c)
                        .template getRealParams<Opm::EclMultiplexerApproach::Default>();
                appendPoints(key, dp.gasOilParams().imbibitionParams().scaledPoints());
                appendPoints(key, dp.oilWaterParams().imbibitionParams().scaledPoints());
                // The kr vertical factors below are derived from the drainage
                // endpoints and would be applied to an imbibition value too,
                // so key them instead of factoring them. pc stays factored —
                // it is not hysteretic, which is what keeps SWATINIT cheap.
                const auto qk = [](double v) {
                    return static_cast<std::int64_t>(std::llround(v * 1.0e9));
                };
                for (double v : { info.maxKrw, info.maxKrow, info.maxKrg, info.maxKrog }) {
                    key.push_back(qk(v));
                }
            }
            // Diagnostic: one table per cell, to tell a dedup/vertical-factor
            // bug apart from a wrong evaluation.
            if (std::getenv("OPM_BAKED_NO_DEDUP") != nullptr) {
                key.push_back(static_cast<std::int64_t>(c));
            }
            const std::array<double, numFactors> refs =
                { info.maxKrw, info.maxKrow, info.maxPcow,
                  info.maxKrg, info.maxKrog, info.maxPcgo };
            auto [it, inserted] = dedup.try_emplace(key, static_cast<std::uint32_t>(tables_.size()));
            if (inserted) {
                tables_.push_back(bakeOne(manager, c, info.Swl));
                tableRefs.push_back(refs);
            }
            cellTable_[c] = it->second;
            cellSwco_[c] = tables_[it->second].swco;
            for (int k = 0; k < numFactors; ++k) {
                cellFactors_[numFactors * c + k] = linFactor(refs[k], tableRefs[it->second][k]);
            }
        }

        if constexpr (hasHysteresis) {
            refreshHysteresisState(manager, nc);
        }
        spotValidate(manager, nc);
    }

    // Reads the per-cell Carlson state from the manager. The tables never
    // change; only these four scalars per cell do, so this is what a
    // simulator would call once per accepted time step.
    void refreshHysteresisState(const Manager& manager, std::size_t nc)
    {
        if constexpr (hasHysteresis) {
            for (std::size_t c = 0; c < nc; ++c) {
                const auto& mp = manager.materialLawParams(c);
                const DefaultParams& dp =
                    mp.template getRealParams<Opm::EclMultiplexerApproach::Default>();
                const double swco = cellSwco_[c];
                // gas-oil: law argument is 1 - Swl - sg, so the condition
                // arg <= krnSwMdc becomes sg >= 1 - Swl - krnSwMdc, and the
                // law's (arg + delta) is (sg - delta) in the sg coordinate.
                cellHyst_[4 * c + 0] = 1.0 - swco - dp.gasOilParams().krnSwMdc();
                cellHyst_[4 * c + 1] = dp.gasOilParams().deltaSwImbKrn();
                // oil-water: law argument is S itself
                cellHyst_[4 * c + 2] = dp.oilWaterParams().krnSwMdc();
                cellHyst_[4 * c + 3] = dp.oilWaterParams().deltaSwImbKrn();
            }
        }
    }

    // Search-free variant (uniform re-tabulation); same interface as evaluate().
    template <class Evaluation>
    void evaluateUniform(std::size_t cellIdx,
                         const Evaluation& sw,
                         const Evaluation& sg,
                         std::array<Evaluation, 3>& kr,
                         Evaluation& pcowDiff,
                         Evaluation& pcgoDiff) const
    {
        const BakedTable& t = tables_[cellTable_[cellIdx]];
        const double swco = cellSwco_[cellIdx];
        const double* f = &cellFactors_[numFactors * cellIdx];

        Evaluation krw, krg, kroOw, kroGo;
        t.uSw.eval2(sw, krw, pcowDiff);
        krw *= f[0];
        pcowDiff *= f[2];
        t.uSg.eval2(sg, krg, pcgoDiff);
        krg *= f[3];
        pcgoDiff *= f[5];

        const Evaluation swClamped = max(Evaluation(swco), sw);
        const Evaluation S = sg + swClamped;
        t.uSo.eval2(S, kroOw, kroGo);
        kroOw *= f[1];
        kroGo *= f[4];

        Evaluation kro;
        constexpr double epsilon = 1e-5;
        const double Sval = Opm::getValue(S);
        if (Sval - swco < epsilon) {
            const Evaluation kro2 = (kroOw + kroGo) / 2;
            if (Sval - swco > epsilon / 2) {
                const Evaluation kro1 = (sg * kroGo + (swClamped - swco) * kroOw) / (S - swco);
                const Evaluation alpha = (epsilon - (S - swco)) / (epsilon / 2);
                kro = kro2 * alpha + kro1 * (1 - alpha);
            }
            else {
                kro = kro2;
            }
        }
        else {
            kro = (sg * kroGo + (swClamped - swco) * kroOw) / (S - swco);
        }
        kr[waterIdx] = krw;
        kr[oilIdx] = kro;
        kr[gasIdx] = krg;
    }

    // Benchmark aid only: drive every cell onto the imbibition branch so the
    // reversed path can be priced. Not a physical state.
    void forceReversedForBenchmark()
    {
        if constexpr (hasHysteresis) {
            for (std::size_t c = 0; c * 4 < cellHyst_.size(); ++c) {
                cellHyst_[4 * c + 0] = -1.0;   // sg  > thr always
                cellHyst_[4 * c + 1] = 0.0;
                cellHyst_[4 * c + 2] = -1.0;   // S   > thr always
                cellHyst_[4 * c + 3] = 0.0;
            }
        }
    }

    std::size_t numTables() const { return tables_.size(); }

    // Storage split into the randomly-accessed table pool (must stay
    // cache-resident) and the sequentially-streamed per-cell arrays.
    struct Storage
    {
        std::size_t poolBytes = 0;      // node arrays, random access
        std::size_t perCellBytes = 0;   // index + swco + factors, streamed
        std::size_t distinctFactorSets = 0;
        std::size_t nc = 0;
    };

    Storage storage() const
    {
        Storage s;
        s.nc = cellTable_.size();
        for (const auto& t : tables_) {
            for (const BakedCurve* c : { &t.sw, &t.sg, &t.so, &t.sgImb, &t.soImb }) {
                s.poolBytes += c->x.size() * sizeof(double)
                             + c->vals.size() * sizeof(double)
                             + sizeof(c->hint);
            }
        }
        s.perCellBytes = cellTable_.size() * sizeof(std::uint32_t)
                       + cellSwco_.size() * sizeof(double)
                       + cellFactors_.size() * sizeof(double)
                       + cellHyst_.size() * sizeof(double);
        std::set<std::vector<std::int64_t>> sets;
        for (std::size_t c = 0; c < s.nc; ++c) {
            std::vector<std::int64_t> k(numFactors);
            for (int j = 0; j < numFactors; ++j) {
                k[j] = static_cast<std::int64_t>(
                    std::llround(cellFactors_[numFactors * c + j] * 1.0e9));
            }
            sets.insert(std::move(k));
        }
        s.distinctFactorSets = sets.size();
        return s;
    }
    double maxUniformError() const { return maxUniformError_; }
    std::size_t numUniformNodes() const
    {
        std::size_t n = 0;
        for (const auto& t : tables_) {
            n += t.uSw.n + t.uSg.n + t.uSo.n;
        }
        return n;
    }
    std::size_t numNodes() const
    {
        std::size_t n = 0;
        for (const auto& t : tables_) {
            n += t.sw.x.size() + t.sg.x.size() + t.so.x.size();
        }
        return n;
    }
    double maxBakeError() const { return maxBakeError_; }

    // Fused evaluation: kr[3], pcow-diff, pcgo-diff. One search per axis,
    // analytic derivatives, default (saturation-weighted) 3-phase oil combine (with the same
    // epsilon regularization as EclDefaultMaterial::krn).
    template <class Evaluation>
    void evaluate(std::size_t cellIdx,
                  const Evaluation& sw,
                  const Evaluation& sg,
                  std::array<Evaluation, 3>& kr,     // [water, oil, gas] canonical idx below
                  Evaluation& pcowDiff,
                  Evaluation& pcgoDiff) const
    {
        const BakedTable& t = tables_[cellTable_[cellIdx]];
        const double swco = cellSwco_[cellIdx];
        const double* f = &cellFactors_[numFactors * cellIdx];

        Evaluation krw, krg, kroOw, kroGo;
        t.sw.eval2(sw, krw, pcowDiff);
        krw *= f[0];
        pcowDiff *= f[2];
        t.sg.eval2(sg, krg, pcgoDiff);
        pcgoDiff *= f[5];

        const Evaluation swClamped = max(Evaluation(swco), sw);
        const Evaluation S = sg + swClamped;
        t.so.eval2(S, kroOw, kroGo);
        kroGo *= f[4];

        // Carlson: the two non-wetting relperms switch to the imbibition
        // curve, shifted, once the cell has reversed. One compare and one
        // add each; unreversed cells keep the three-lookup path.
        if constexpr (hasHysteresis) {
            if (t.hasImb) {
                const double* h = &cellHyst_[4 * cellIdx];
                Evaluation dummy;
                // gas-oil: the law's argument is 1-Swl-sg, so its
                // "arg <= krnSwMdc means drainage" maps to sg >= threshold;
                // the imbibition branch is therefore sg < threshold.
                // The Carlson shift can push the argument outside the curve's
                // range, where the true law clamps to the end value with zero
                // derivative rather than extrapolating along the end slope.
                if (Opm::getValue(sg) < h[0]) {
                    t.sgImb.eval2(clampToRange(sg - h[1], t.sgImb), krg, dummy);
                }
                // oil-water: the argument is S itself, so the sense is direct.
                if (Opm::getValue(S) > h[2]) {
                    t.soImb.eval2(clampToRange(S + h[3], t.soImb), kroOw, dummy);
                }
            }
        }
        krg *= f[3];
        kroOw *= f[1];   // factors applied before the three-phase combine

        Evaluation kro;
        constexpr double epsilon = 1e-5;
        const double Sval = Opm::getValue(S);
        if (Sval - swco < epsilon) {
            const Evaluation kro2 = (kroOw + kroGo) / 2;
            if (Sval - swco > epsilon / 2) {
                const Evaluation kro1 = (sg * kroGo + (swClamped - swco) * kroOw) / (S - swco);
                const Evaluation alpha = (epsilon - (S - swco)) / (epsilon / 2);
                kro = kro2 * alpha + kro1 * (1 - alpha);
            }
            else {
                kro = kro2;
            }
        }
        else {
            kro = (sg * kroGo + (swClamped - swco) * kroOw) / (S - swco);
        }

        kr[waterIdx] = krw;
        kr[oilIdx] = kro;
        kr[gasIdx] = krg;
    }

private:
    static double linFactor(double cellV, double tableV)
    {
        if (tableV > 0.0 && cellV > 0.0) {
            const double fac = cellV / tableV;
            if (std::isfinite(fac)) {
                return fac;
            }
        }
        return 1.0;
    }

    // Constant outside the curve's range, so the derivative vanishes there —
    // PiecewiseLinearTwoPhaseMaterial returns the end value, it does not
    // extrapolate.
    template <class Evaluation>
    static Evaluation clampToRange(const Evaluation& x, const BakedCurve& c)
    {
        const double v = Opm::getValue(x);
        if (v < c.x.front()) { return Evaluation(c.x.front()); }
        if (v > c.x.back())  { return Evaluation(c.x.back()); }
        return x;
    }

    // Every scaling quantity of one EPS point set, quantized into a key. Used
    // for the imbibition curves, whose vertical scales are not factored out.
    template <class Points>
    static void appendPoints(std::vector<std::int64_t>& key, const Points& p)
    {
        const auto q = [](double v) { return static_cast<std::int64_t>(std::llround(v * 1.0e9)); };
        for (double v : { p.maxPcnw(), p.maxKrw(), p.krwr(), p.maxKrn(), p.krnr() }) {
            key.push_back(q(v));
        }
        for (const auto* arr : { &p.saturationPcPoints(), &p.saturationKrwPoints(),
                                 &p.saturationKrnPoints() }) {
            for (double v : *arr) {
                key.push_back(q(v));
            }
        }
    }

    template <class Info, class EpsConfig>
    static std::vector<std::int64_t> makeKey(int satnum, const Info& i,
                                             const EpsConfig& owCfg, const EpsConfig& goCfg)
    {
        const auto q = [](double v) { return static_cast<std::int64_t>(std::llround(v * 1.0e9)); };
        std::vector<std::int64_t> key
            { satnum,
              q(i.Swl), q(i.Swcr), q(i.Swu), q(i.Sgl), q(i.Sgcr), q(i.Sgu),
              q(i.Sowcr), q(i.Sogcr) };
        // Vertical scales enter only as shape information: with 3pt vertical
        // scaling active the residual/max ratio shapes the curve; with plain
        // 2pt scaling the magnitude is a pure per-cell factor and does not
        // key. A non-positive max cannot be factored and keys the raw pair.
        const auto shape = [&](double r, double m, bool threePt) {
            if (m > 0.0) {
                key.push_back(threePt ? q(r / m) : 0);
            }
            else {
                key.push_back(std::numeric_limits<std::int64_t>::min());
                key.push_back(q(r));
                key.push_back(q(m));
            }
        };
        // gas-oil system: wetting = oil (kro_go), non-wetting = gas (krg)
        shape(i.Krwr, i.maxKrw, owCfg.enableThreePointKrwScaling());
        shape(i.Krorw, i.maxKrow, owCfg.enableThreePointKrnScaling());
        shape(0.0, i.maxPcow, false);
        shape(i.Krgr, i.maxKrg, goCfg.enableThreePointKrnScaling());
        shape(i.Krorg, i.maxKrog, goCfg.enableThreePointKrwScaling());
        shape(0.0, i.maxPcgo, false);
        return key;
    }

    // Compare the baked evaluation (incl. the factored per-cell pc scale)
    // against the true material law on a sample of cells and states.
    void spotValidate(const Manager& manager, std::size_t nc)
    {
        const std::size_t stride = std::max<std::size_t>(nc / 200, 1);
        for (std::size_t c = 0; c < nc; c += stride) {
            const auto& mp = manager.materialLawParams(c);
            const DefaultParams& dp =
                mp.template getRealParams<Opm::EclMultiplexerApproach::Default>();
            for (int a = 0; a <= 3; ++a) {
                for (int b = 0; a + b <= 3; ++b) {
                    const double sw = 0.05 + 0.3 * a;
                    const double sg = 0.02 + 0.3 * b;
                    FS fs;
                    fs.sat[waterIdx] = sw;
                    fs.sat[gasIdx] = sg;
                    fs.sat[oilIdx] = 1.0 - sw - sg;
                    std::array<double, numPhases> krT{}, pcT{};
                    DefaultMaterial::relativePermeabilities(krT, dp, fs);
                    DefaultMaterial::capillaryPressures(pcT, dp, fs);
                    std::array<double, 3> krB;
                    double pcow, pcgo;
                    evaluate(c, sw, sg, krB, pcow, pcgo);
                    double err = 0.0;
                    for (int ph = 0; ph < numPhases; ++ph) {
                        err = std::max(err, std::abs(krB[ph] - krT[ph])
                                            / std::max(std::abs(krT[ph]), 1.0));
                    }
                    const double pcowT = pcT[waterIdx] - pcT[oilIdx];
                    const double pcgoT = pcT[gasIdx] - pcT[oilIdx];
                    err = std::max(err, std::abs(pcow - pcowT) / std::max(std::abs(pcowT), 1.0));
                    err = std::max(err, std::abs(pcgo - pcgoT) / std::max(std::abs(pcgoT), 1.0));
                    maxBakeError_ = std::max(maxBakeError_, err);
                }
            }
        }
    }

    template <int NC = 2, class Fn>
    BakedCurve bakeCurve(const Fn& truth, std::set<double>& nodes)
    {
        // refine until every interval's midpoint matches the true law
        constexpr double tol = 1.0e-12;
        for (int pass = 0; pass < 12; ++pass) {
            std::vector<double> xs(nodes.begin(), nodes.end());
            std::vector<double> insert;
            for (std::size_t i = 0; i + 1 < xs.size(); ++i) {
                const double xm = 0.5 * (xs[i] + xs[i + 1]);
                if (xs[i + 1] - xs[i] < 1.0e-10) {
                    continue;
                }
                const auto lo = truth(xs[i]);
                const auto hi = truth(xs[i + 1]);
                const auto mid = truth(xm);
                for (int col = 0; col < NC; ++col) {
                    const double lin = 0.5 * (lo[col] + hi[col]);
                    const double scale = std::max({ std::abs(lo[col]), std::abs(hi[col]), 1.0 });
                    if (std::abs(mid[col] - lin) > tol * scale) {
                        insert.push_back(xm);
                        // kink snap: for piecewise-linear truth with one kink in
                        // the interval, intersecting the secants through the
                        // left and right quarters yields the kink exactly.
                        const double w = xs[i + 1] - xs[i];
                        const double xq1 = xs[i] + 0.25 * w;
                        const double xq3 = xs[i] + 0.75 * w;
                        const auto q1 = truth(xq1);
                        const auto q3 = truth(xq3);
                        const double sl = (q1[col] - lo[col]) / (xq1 - xs[i]);
                        const double sr = (hi[col] - q3[col]) / (xs[i + 1] - xq3);
                        if (std::abs(sl - sr) > 1.0e-14 * std::max(std::abs(sl), std::abs(sr))) {
                            // lines: lo + sl (x - xi)  and  hi + sr (x - xi1)
                            const double xk = (hi[col] - lo[col] + sl * xs[i] - sr * xs[i + 1])
                                              / (sl - sr);
                            if (xk > xs[i] + 1.0e-12 && xk < xs[i + 1] - 1.0e-12) {
                                insert.push_back(xk);
                            }
                        }
                        break;
                    }
                }
            }
            if (insert.empty()) {
                break;
            }
            nodes.insert(insert.begin(), insert.end());
        }

        // sample final nodes, then prune collinear interior nodes
        std::vector<double> xs(nodes.begin(), nodes.end());
        std::vector<std::array<double, NC>> ys(xs.size());
        for (std::size_t i = 0; i < xs.size(); ++i) {
            ys[i] = truth(xs[i]);
        }
        std::vector<char> keep(xs.size(), 1);
        for (std::size_t i = 1; i + 1 < xs.size(); ++i) {
            bool collinear = true;
            // find previous kept node
            std::size_t p = i - 1;
            while (!keep[p]) { --p; }
            const double t = (xs[i] - xs[p]) / (xs[i + 1] - xs[p]);
            for (int col = 0; col < NC && collinear; ++col) {
                const double lin = ys[p][col] * (1.0 - t) + ys[i + 1][col] * t;
                const double scale = std::max(std::abs(ys[i][col]), 1.0);
                collinear = std::abs(lin - ys[i][col]) <= 1.0e-13 * scale;
            }
            if (collinear) {
                keep[i] = 0;
            }
        }

        BakedCurve c;
        c.nCols = NC;
        std::vector<std::array<double, NC>> kv;
        for (std::size_t i = 0; i < xs.size(); ++i) {
            if (keep[i]) {
                c.x.push_back(xs[i]);
                kv.push_back(ys[i]);
            }
        }
        const std::size_t n = c.x.size();
        c.vals.assign(2 * NC * n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const double idx = (i + 1 < n) ? 1.0 / (c.x[i + 1] - c.x[i]) : 0.0;
            for (int col = 0; col < NC; ++col) {
                c.vals[2 * NC * i + 2 * col] = kv[i][col];
                c.vals[2 * NC * i + 2 * col + 1] =
                    (i + 1 < n) ? (kv[i + 1][col] - kv[i][col]) * idx : 0.0;
            }
        }

        // validation probe (golden-ratio points; belt and braces)
        double maxErr = 0.0;
        for (std::size_t i = 0; i + 1 < n; ++i) {
            const double xp = c.x[i] + 0.618033988749895 * (c.x[i + 1] - c.x[i]);
            const auto tv = truth(xp);
            const double dx = xp - c.x[i];
            for (int col = 0; col < NC; ++col) {
                const double pv = c.vals[2 * NC * i + 2 * col]
                                + c.vals[2 * NC * i + 2 * col + 1] * dx;
                maxErr = std::max(maxErr,
                                  std::abs(pv - tv[col]) / std::max(std::abs(tv[col]), 1.0));
            }
        }
        maxBakeError_ = std::max(maxBakeError_, maxErr);
        c.buildHint();
        return c;
    }

    // Re-tabulate on a uniform abscissa, doubling the node count until the
    // deviation from the true curve is under uniformTol_ (relative).
    template <class Fn>
    UniformCurve bakeUniform(const Fn& truth)
    {
        UniformCurve u;
        for (std::size_t n = 65; n <= 8193; n = 2 * (n - 1) + 1) {
            u.n = n;
            u.x0 = 0.0;
            u.invH = double(n - 1);
            u.vals.assign(4 * n, 0.0);
            std::vector<std::array<double, 2>> y(n);
            for (std::size_t i = 0; i < n; ++i) {
                y[i] = truth(double(i) / u.invH);
            }
            for (std::size_t i = 0; i < n; ++i) {
                const bool last = (i + 1 == n);
                u.vals[4 * i + 0] = y[i][0];
                u.vals[4 * i + 1] = last ? 0.0 : (y[i + 1][0] - y[i][0]) * u.invH;
                u.vals[4 * i + 2] = y[i][1];
                u.vals[4 * i + 3] = last ? 0.0 : (y[i + 1][1] - y[i][1]) * u.invH;
            }
            // probe interval interiors for the worst deviation
            double err = 0.0;
            for (std::size_t i = 0; i + 1 < n; ++i) {
                for (const double frac : { 0.3183098861837907, 0.618033988749895 }) {
                    const double x = (double(i) + frac) / u.invH;
                    const auto tv = truth(x);
                    const double dx = x - double(i) / u.invH;
                    const double p0 = u.vals[4 * i + 0] + u.vals[4 * i + 1] * dx;
                    const double p1 = u.vals[4 * i + 2] + u.vals[4 * i + 3] * dx;
                    err = std::max(err, std::abs(p0 - tv[0]) / std::max(std::abs(tv[0]), 1.0));
                    err = std::max(err, std::abs(p1 - tv[1]) / std::max(std::abs(tv[1]), 1.0));
                }
            }
            maxUniformError_ = std::max(maxUniformError_, err);
            if (err <= uniformTol_) {
                break;
            }
        }
        return u;
    }

    BakedTable bakeOne(const Manager& manager, std::size_t cellIdx, double swco)
    {
        const auto& mp = manager.materialLawParams(cellIdx);
        const auto& info = manager.oilWaterScaledEpsInfoDrainage(cellIdx);
        const auto seed = [&](std::initializer_list<double> anchors) {
            std::set<double> s;
            for (int i = 0; i <= 64; ++i) {
                s.insert(double(i) / 64.0);
            }
            for (double a : anchors) {
                if (a > 0.0 && a < 1.0) {
                    s.insert(a);
                }
            }
            return s;
        };

        const DefaultParams& dp =
            mp.template getRealParams<Opm::EclMultiplexerApproach::Default>();
        Truth truth{ &dp, swco };

        // candidate nodes: endpoints + scaled EPS anchors + uniform fill;
        // the refine loop discovers any interior table kinks.

        BakedTable t;
        t.swco = swco;
        {
            auto nodes = seed({ info.Swl, info.Swcr, info.Swu, 1.0 - info.Sowcr,
                                1.0 - info.Sowcr - info.Sgl });
            const auto f = [&](double x) { return truth.swAxis(x); };
            t.sw = bakeCurve(f, nodes);
            t.uSw = bakeUniform(f);
        }
        {
            auto nodes = seed({ info.Sgl, info.Sgcr, info.Sgu, 1.0 - info.Sogcr - info.Swl });
            const auto f = [&](double x) { return truth.sgAxis(x); };
            t.sg = bakeCurve(f, nodes);
            t.uSg = bakeUniform(f);
        }
        if constexpr (hasHysteresis) {
            // Seeded from the *imbibition* endpoints (ISGCR & co): its kinks
            // sit elsewhere than the drainage ones, and bake-refine only
            // guarantees the tolerance at seeded kinks.
            const auto& gip = truth.p->gasOilParams().imbibitionParams().scaledPoints();
            auto nodes = seed({ info.Sgl, info.Sgcr, info.Sgu, 1.0 - info.Sogcr - info.Swl,
                                1.0 - info.Swl - gip.saturationKrnPoints()[0],
                                1.0 - info.Swl - gip.saturationKrnPoints()[1],
                                1.0 - info.Swl - gip.saturationKrnPoints()[2] });
            t.sgImb = bakeCurve([&](double x) { return truth.sgAxisImb(x); }, nodes);
            t.hasImb = true;
        }
        {
            auto nodes = seed({ info.Swl, info.Swcr, 1.0 - info.Sowcr, info.Swl + info.Sgcr,
                                1.0 - info.Sogcr });
            const auto f = [&](double x) { return truth.soAxis(x); };
            t.so = bakeCurve(f, nodes);
            t.uSo = bakeUniform(f);
        }
        if constexpr (hasHysteresis) {
            const auto& oip = truth.p->oilWaterParams().imbibitionParams().scaledPoints();
            auto nodes = seed({ info.Swl, info.Swcr, 1.0 - info.Sowcr, info.Swl + info.Sgcr,
                                1.0 - info.Sogcr,
                                oip.saturationKrnPoints()[0], oip.saturationKrnPoints()[1],
                                oip.saturationKrnPoints()[2] });
            t.soImb = bakeCurve([&](double x) { return truth.soAxisImb(x); }, nodes);
        }
        return t;
    }

    std::vector<BakedTable> tables_;
    std::vector<std::uint32_t> cellTable_;
    std::vector<double> cellSwco_;
    std::vector<double> cellFactors_;   // numFactors per cell
    // Carlson state, 2 per system: threshold in baked coordinate + shift.
    // krg  : drainage while sg <= sgThr,  else imbibition curve at sg - sgShift
    // kroOw: drainage while S  <= sThr,   else imbibition curve at S  + sShift
    std::vector<double> cellHyst_;      // 4 per cell: sgThr, sgShift, sThr, sShift
    std::vector<double> cellStone_;     // 2 per cell: krocw, eta (Stone only)
    double maxBakeError_ = 0.0;
    double maxUniformError_ = 0.0;
    double uniformTol_ = 1.0e-6;
};

} // namespace Opm::Bench

#endif
