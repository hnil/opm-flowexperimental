// Standalone intensive-quantities benchmark on a real deck, no simulator.
//
// Compares, on identical fabricated states:
//   A  master BlackOilIntensiveQuantities::update(problem, priVars, i, 0)
//   B  BlackOilIntensiveQuantitiesFast::update(...)          (flowexperimental)
//   C  decomposed region-batched passes (SATNUM relperm, PVTNUM pvt with
//      concrete hoisted PVT objects), plus generic sub-part loops that give
//      the PVT/relperm split and the cell-order (function-switching) cost.
//
// Usage: bench_iq_evaluation DECK [--n N] [--undersat-frac F] [--seed S] [--reps R]
#include "config.h"

#include "../../../opm-simulators/tests/TestTypeTag.hpp"

#include <opm/models/blackoil/blackoilintensivequantities.hh>
#include <opm/flowexperimental/blackoilintensivequantitiesfast.hh>

#include "StubIQProblem.hpp"
#include "BakedRelperm.hpp"

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidmatrixinteractions/EclDefaultMaterial.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/InputErrorAction.hpp>
#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Python/Python.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

// Stub problem variant providing the baked-relperm hook, so the fast IQ's
// compile-time detection activates: full IQ update, baked relperm inside the
// same mixed per-cell loop (the simulator situation).
template <class TypeTag>
class StubIQProblemBaked : public Opm::Bench::StubIQProblem<TypeTag>
{
    using Base = Opm::Bench::StubIQProblem<TypeTag>;
    using MaterialLawT = Opm::GetPropType<TypeTag, Opm::Properties::MaterialLaw>;
    using ScalarT = Opm::GetPropType<TypeTag, Opm::Properties::Scalar>;
    using FluidSystemT = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using ManagerT = Opm::EclMaterialLaw::Manager<typename MaterialLawT::Traits>;

public:
    using BakedT = Opm::Bench::BakedRelperm<MaterialLawT, ManagerT, ScalarT>;
    using Base::Base;

    void setBaked(const BakedT* b) { baked_ = b; }

    template <class Ev>
    void evalBakedRelperm(unsigned globalSpaceIdx,
                          const Ev& sw,
                          const Ev& sg,
                          std::array<Ev, FluidSystemT::numPhases>& mobility,
                          std::array<Ev, FluidSystemT::numPhases>& pC) const
    {
        baked_->evaluate(globalSpaceIdx, sw, sg, mobility,
                         pC[FluidSystemT::waterPhaseIdx], pC[FluidSystemT::gasPhaseIdx]);
        pC[FluidSystemT::oilPhaseIdx] = 0.0;
    }

private:
    const BakedT* baked_ = nullptr;
};

// ---------------------------------------------------------------------------
// TypeTags: identical except the IntensiveQuantities implementation.
// ---------------------------------------------------------------------------
namespace Opm::Properties {
namespace TTag {
struct BenchIQMaster { using InheritsFrom = std::tuple<TestTypeTag>; };
struct BenchIQFast   { using InheritsFrom = std::tuple<TestTypeTag>; };
struct BenchIQFastBaked { using InheritsFrom = std::tuple<TestTypeTag>; };
} // namespace TTag

template <class TypeTag, class MyTypeTag>
struct BenchMaterialLaw
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Traits = ThreePhaseMaterialTraits<Scalar,
                                            FluidSystem::waterPhaseIdx,
                                            FluidSystem::oilPhaseIdx,
                                            FluidSystem::gasPhaseIdx,
                                            /*enableHysteresis=*/false,
                                            /*enableEndpointScaling=*/true>;
public:
    using EclMaterialLawManager = ::Opm::EclMaterialLaw::Manager<Traits>;
    using type = typename EclMaterialLawManager::MaterialLaw;
};

template <class TypeTag>
struct Problem<TypeTag, TTag::BenchIQMaster> { using type = ::Opm::Bench::StubIQProblem<TypeTag>; };
template <class TypeTag>
struct Problem<TypeTag, TTag::BenchIQFast> { using type = ::Opm::Bench::StubIQProblem<TypeTag>; };
template <class TypeTag>
struct Problem<TypeTag, TTag::BenchIQFastBaked> { using type = ::StubIQProblemBaked<TypeTag>; };

template <class TypeTag>
struct MaterialLaw<TypeTag, TTag::BenchIQMaster> : BenchMaterialLaw<TypeTag, TTag::BenchIQMaster> {};
template <class TypeTag>
struct MaterialLaw<TypeTag, TTag::BenchIQFast> : BenchMaterialLaw<TypeTag, TTag::BenchIQFast> {};
template <class TypeTag>
struct MaterialLaw<TypeTag, TTag::BenchIQFastBaked> : BenchMaterialLaw<TypeTag, TTag::BenchIQFastBaked> {};

// The context-free updates assert these off.
template <class TypeTag>
struct EnableDiffusion<TypeTag, TTag::BenchIQMaster> { static constexpr bool value = false; };
template <class TypeTag>
struct EnableDiffusion<TypeTag, TTag::BenchIQFast> { static constexpr bool value = false; };
template <class TypeTag>
struct EnableDiffusion<TypeTag, TTag::BenchIQFastBaked> { static constexpr bool value = false; };
template <class TypeTag>
struct EnableDispersion<TypeTag, TTag::BenchIQMaster> { static constexpr bool value = false; };
template <class TypeTag>
struct EnableDispersion<TypeTag, TTag::BenchIQFast> { static constexpr bool value = false; };
template <class TypeTag>
struct EnableDispersion<TypeTag, TTag::BenchIQFastBaked> { static constexpr bool value = false; };

template <class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::BenchIQFast>
{ using type = BlackOilIntensiveQuantitiesFast<TypeTag>; };
template <class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::BenchIQFastBaked>
{ using type = BlackOilIntensiveQuantitiesFast<TypeTag>; };

} // namespace Opm::Properties

namespace {

using MasterTag = Opm::Properties::TTag::BenchIQMaster;
using FastTag = Opm::Properties::TTag::BenchIQFast;

using Scalar = Opm::GetPropType<MasterTag, Opm::Properties::Scalar>;
using FluidSystem = Opm::GetPropType<MasterTag, Opm::Properties::FluidSystem>;
using MaterialLaw = Opm::GetPropType<MasterTag, Opm::Properties::MaterialLaw>;
using Manager = Opm::EclMaterialLaw::Manager<typename MaterialLaw::Traits>;

constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;
constexpr unsigned oilPhaseIdx = FluidSystem::oilPhaseIdx;
constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
constexpr int numPhases = FluidSystem::numPhases;

constexpr Scalar benchTemperature = 323.15;

// ---------------------------------------------------------------------------
// CLI + small helpers
// ---------------------------------------------------------------------------
struct Options
{
    std::string deck;
    std::size_t n = 0;          // 0 = all active cells
    double undersatFrac = 0.3;
    unsigned seed = 12345;
    int reps = 10;
    bool checkOnly = false;     // report baked-path coverage and exit
};

Options parseArgs(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: bench_iq_evaluation DECK [--n N] [--undersat-frac F]"
                     " [--seed S] [--reps R]\n";
        std::exit(1);
    }
    Options o;
    o.deck = argv[1];
    for (int i = 2; i + 1 < argc; i += 2) {
        const std::string k = argv[i];
        const std::string v = argv[i + 1];
        if (k == "--n") { o.n = std::stoul(v); }
        else if (k == "--undersat-frac") { o.undersatFrac = std::stod(v); }
        else if (k == "--seed") { o.seed = static_cast<unsigned>(std::stoul(v)); }
        else if (k == "--reps") { o.reps = std::stoi(v); }
        else if (k == "--check-only") { o.checkOnly = (v != "0"); }
        else { std::cerr << "unknown option " << k << "\n"; std::exit(1); }
    }
    return o;
}

double toNsPerCell(std::chrono::steady_clock::duration d, std::size_t nc)
{
    return std::chrono::duration<double, std::nano>(d).count() / double(nc);
}

// Fabricated per-cell state; consumed by every variant.
struct FabState
{
    Scalar po;
    Scalar sw;
    Scalar sg;      // only meaningful when !undersat
    Scalar rs;      // only meaningful when undersat
    bool undersat;
};

std::vector<FabState> fabricateStates(std::size_t nc,
                                      const std::vector<int>& pvtnum,
                                      double undersatFrac,
                                      unsigned seed)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> upo(150.0e5, 320.0e5);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    std::uniform_real_distribution<double> ufrac(0.3, 0.95);

    std::vector<FabState> states(nc);
    for (std::size_t i = 0; i < nc; ++i) {
        FabState& s = states[i];
        s.po = upo(gen);
        s.undersat = FluidSystem::enableDissolvedGas() && (u01(gen) < undersatFrac);
        if (s.undersat) {
            // no free gas: Sw + So = 1, Rs below the saturated value
            s.sw = 0.7 * u01(gen);
            s.sg = 0.0;
            const Scalar rsSat =
                FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtnum[i],
                                                                    Scalar(benchTemperature),
                                                                    s.po);
            s.rs = ufrac(gen) * rsSat;
        }
        else {
            double a = u01(gen);
            double b = u01(gen);
            if (a + b > 1.0) { a = 1.0 - a; b = 1.0 - b; }
            s.sw = a;
            s.sg = b;
            s.rs = 0.0;
        }
    }
    return states;
}

template <class TypeTag>
std::vector<Opm::GetPropType<TypeTag, Opm::Properties::PrimaryVariables>>
makePrimaryVariables(const std::vector<FabState>& states, const std::vector<int>& pvtnum)
{
    using PV = Opm::GetPropType<TypeTag, Opm::Properties::PrimaryVariables>;
    using Indices = Opm::GetPropType<TypeTag, Opm::Properties::Indices>;

    std::vector<PV> pvs(states.size());
    for (std::size_t i = 0; i < states.size(); ++i) {
        PV& pv = pvs[i];
        for (std::size_t k = 0; k < pv.size(); ++k) {
            pv[k] = 0.0;
        }
        pv.setPrimaryVarsMeaningWater(PV::WaterMeaning::Sw);
        pv.setPrimaryVarsMeaningPressure(PV::PressureMeaning::Po);
        pv.setPrimaryVarsMeaningBrine(PV::BrineMeaning::Disabled);
        pv.setPrimaryVarsMeaningSolvent(PV::SolventMeaning::Disabled);
        pv.setPvtRegionIndex(pvtnum[i]);
        pv[Indices::pressureSwitchIdx] = states[i].po;
        pv[Indices::waterSwitchIdx] = states[i].sw;
        if (states[i].undersat) {
            pv.setPrimaryVarsMeaningGas(PV::GasMeaning::Rs);
            pv[Indices::compositionSwitchIdx] = states[i].rs;
        }
        else {
            pv.setPrimaryVarsMeaningGas(PV::GasMeaning::Sg);
            pv[Indices::compositionSwitchIdx] = states[i].sg;
        }
    }
    return pvs;
}

// Checksum over the quantities every variant can produce: phase mobilities,
// inverse formation volume factors and capillary-pressure differences.
template <class IQ>
long double iqChecksum(const IQ& iq)
{
    long double c = 0.0;
    for (int ph = 0; ph < numPhases; ++ph) {
        if (!FluidSystem::phaseIsActive(ph)) {
            continue;
        }
        c += Opm::getValue(iq.mobility(ph));
        c += Opm::getValue(iq.fluidState().invB(ph));
    }
    const auto& fs = iq.fluidState();
    c += Opm::getValue(fs.pressure(gasPhaseIdx)) - Opm::getValue(fs.pressure(oilPhaseIdx));
    c += Opm::getValue(fs.pressure(waterPhaseIdx)) - Opm::getValue(fs.pressure(oilPhaseIdx));
    return c;
}

struct RunResult
{
    double nsPerCell;
    long double checksum;
};

template <class TypeTag>
RunResult runIQVariant(const Opm::GetPropType<TypeTag, Opm::Properties::Problem>& problem,
                       const std::vector<Opm::GetPropType<TypeTag, Opm::Properties::PrimaryVariables>>& pvs,
                       const std::vector<std::size_t>& order,
                       int reps)
{
    using IQ = Opm::GetPropType<TypeTag, Opm::Properties::IntensiveQuantities>;

    volatile double sink = 0.0;
    long double checksum = 0.0;
    auto best = std::chrono::steady_clock::duration::max();
    for (int rep = 0; rep < reps; ++rep) {
        long double chk = 0.0;
        const auto t0 = std::chrono::steady_clock::now();
        for (const std::size_t i : order) {
            IQ iq;
            iq.update(problem, pvs[i], static_cast<unsigned>(i), /*timeIdx=*/0);
            chk += iqChecksum(iq);
        }
        const auto dt = std::chrono::steady_clock::now() - t0;
        best = std::min(best, dt);
        sink = sink + static_cast<double>(chk);
        if (rep == 0) {
            checksum = chk;
        }
    }
    (void)sink;
    return { toNsPerCell(best, order.size()), checksum };
}

// Resolved states (from a master-IQ sweep) feeding the sub-part loops so all
// parts see identical, IQ-consistent inputs.
struct ResolvedCell
{
    Scalar pw, po, pg;
    Scalar sw, so, sg;
    Scalar rs, rv;
};

template <class TypeTag>
std::vector<ResolvedCell>
resolveStates(const Opm::Bench::StubIQProblem<TypeTag>& problem,
              const std::vector<Opm::GetPropType<TypeTag, Opm::Properties::PrimaryVariables>>& pvs)
{
    using IQ = Opm::GetPropType<TypeTag, Opm::Properties::IntensiveQuantities>;

    std::vector<ResolvedCell> res(pvs.size());
    for (std::size_t i = 0; i < pvs.size(); ++i) {
        IQ iq;
        iq.update(problem, pvs[i], static_cast<unsigned>(i), 0);
        const auto& fs = iq.fluidState();
        res[i] = { Opm::getValue(fs.pressure(waterPhaseIdx)),
                   Opm::getValue(fs.pressure(oilPhaseIdx)),
                   Opm::getValue(fs.pressure(gasPhaseIdx)),
                   Opm::getValue(fs.saturation(waterPhaseIdx)),
                   Opm::getValue(fs.saturation(oilPhaseIdx)),
                   Opm::getValue(fs.saturation(gasPhaseIdx)),
                   Opm::getValue(fs.Rs()),
                   Opm::getValue(fs.Rv()) };
    }
    return res;
}

// Build a fluid state for the sub-part loops. For AD types the primary
// unknowns (po, Sw, Sg) carry derivatives, mirroring the simulator.
template <class Evaluation>
Opm::BlackOilFluidState<Evaluation, FluidSystem>
makeFluidState(const ResolvedCell& r, unsigned pvtRegionIdx)
{
    using Toolbox = Opm::MathToolbox<Evaluation>;
    Opm::BlackOilFluidState<Evaluation, FluidSystem> fs;
    fs.setPvtRegionIndex(pvtRegionIdx);

    Evaluation po, sw, sg;
    if constexpr (std::is_same_v<Evaluation, Scalar>) {
        po = r.po; sw = r.sw; sg = r.sg;
    }
    else {
        po = Toolbox::createVariable(r.po, 0);
        sw = Toolbox::createVariable(r.sw, 1);
        sg = Toolbox::createVariable(r.sg, 2);
    }
    const Evaluation so = 1.0 - sw - sg;
    fs.setSaturation(waterPhaseIdx, sw);
    fs.setSaturation(oilPhaseIdx, so);
    fs.setSaturation(gasPhaseIdx, sg);
    fs.setPressure(oilPhaseIdx, po);
    fs.setPressure(waterPhaseIdx, po + (r.pw - r.po));
    fs.setPressure(gasPhaseIdx, po + (r.pg - r.po));
    fs.setRs(r.rs);
    fs.setRv(r.rv);
    return fs;
}

// Pre-built fluid states: in the simulator the IQ's fluid state exists before
// relperm/PVT run, so construction must stay outside the timed loops.
template <class Evaluation>
std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>>
buildFluidStates(const std::vector<ResolvedCell>& cells, const std::vector<int>& pvtnum)
{
    std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>> fss;
    fss.reserve(cells.size());
    for (std::size_t i = 0; i < cells.size(); ++i) {
        fss.push_back(makeFluidState<Evaluation>(cells[i], pvtnum[i]));
    }
    return fss;
}

// Generic relperm-only loop (per-cell multiplexer params, as the IQ does it).
template <class Evaluation, class LawT = MaterialLaw, class MgrT = Manager>
RunResult runRelperm(const MgrT& manager,
                     const std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>>& fss,
                     const std::vector<std::size_t>& order,
                     int reps,
                     bool pinParams = false)
{
    using EvalArr = std::array<Evaluation, numPhases>;

    volatile double sink = 0.0;
    long double checksum = 0.0;
    auto best = std::chrono::steady_clock::duration::max();
    for (int rep = 0; rep < reps; ++rep) {
        long double chk = 0.0;
        const auto t0 = std::chrono::steady_clock::now();
        for (const std::size_t i : order) {
            const auto& fs = fss[i];
            // pinParams: same params object for every cell, so the parameter
            // footprint drops out and only the arithmetic/dispatch remains.
            const auto& mp = manager.materialLawParams(pinParams ? 0 : i);
            EvalArr kr, pC;
            LawT::relativePermeabilities(kr, mp, fs);
            LawT::capillaryPressures(pC, mp, fs);
            for (int ph = 0; ph < numPhases; ++ph) {
                chk += Opm::getValue(kr[ph]);
            }
            chk += Opm::getValue(pC[gasPhaseIdx]) - Opm::getValue(pC[oilPhaseIdx]);
            chk += Opm::getValue(pC[waterPhaseIdx]) - Opm::getValue(pC[oilPhaseIdx]);
        }
        const auto dt = std::chrono::steady_clock::now() - t0;
        best = std::min(best, dt);
        sink = sink + static_cast<double>(chk);
        if (rep == 0) {
            checksum = chk;
        }
    }
    (void)sink;
    return { toNsPerCell(best, order.size()), checksum };
}

// Same relperm computation, per-cell parameter object padded to a chosen
// size. Isolates the cost of the per-cell parameter footprint from the
// arithmetic: sizeof(EclDefaultMaterial::Params) is 304 B without hysteresis
// and 1440 B with it, so this prices what a footprint refactor would buy.
// Master's own layout is one heap allocation per cell behind a shared_ptr<void>
// in the multiplexer, so this also measures what pooling them into one
// contiguous arena is worth, independently of the struct size.
template <class Evaluation, int PadBytes, class LawT = MaterialLaw, class MgrT = Manager>
RunResult runRelpermPadded(const MgrT& manager,
                           const std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>>& fss,
                           const std::vector<std::size_t>& order,
                           int reps)
{
    using DefaultMat = Opm::EclDefaultMaterial<typename LawT::Traits,
                                               typename LawT::GasOilMaterialLaw,
                                               typename LawT::OilWaterMaterialLaw>;
    using P = typename DefaultMat::Params;
    static constexpr int extra = PadBytes > int(sizeof(P)) ? PadBytes - int(sizeof(P)) : 1;
    struct Padded { P p; char pad[extra]; };

    std::vector<Padded> arr(fss.size());
    for (std::size_t i = 0; i < fss.size(); ++i) {
        arr[i].p = manager.materialLawParams(i)
                       .template getRealParams<Opm::EclMultiplexerApproach::Default>();
        std::fill(std::begin(arr[i].pad), std::end(arr[i].pad), char(1));
    }

    using EvalArr = std::array<Evaluation, numPhases>;
    volatile double sink = 0.0;
    long double checksum = 0.0;
    auto best = std::chrono::steady_clock::duration::max();
    for (int rep = 0; rep < reps; ++rep) {
        long double chk = 0.0;
        const auto t0 = std::chrono::steady_clock::now();
        for (const std::size_t i : order) {
            const auto& fs = fss[i];
            EvalArr kr, pC;
            DefaultMat::relativePermeabilities(kr, arr[i].p, fs);
            DefaultMat::capillaryPressures(pC, arr[i].p, fs);
            for (int ph = 0; ph < numPhases; ++ph) {
                chk += Opm::getValue(kr[ph]);
            }
            chk += Opm::getValue(pC[gasPhaseIdx]) - Opm::getValue(pC[oilPhaseIdx]);
            chk += Opm::getValue(pC[waterPhaseIdx]) - Opm::getValue(pC[oilPhaseIdx]);
        }
        const auto dt = std::chrono::steady_clock::now() - t0;
        best = std::min(best, dt);
        sink = sink + static_cast<double>(chk);
        if (rep == 0) { checksum = chk; }
    }
    (void)sink;
    return { toNsPerCell(best, order.size()), checksum };
}

// Generic PVT-only loop (FluidSystem multiplexer path, as master IQ).
template <class Evaluation>
RunResult runPvt(const std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>>& fss,
                 const std::vector<std::size_t>& order,
                 int reps)
{
    volatile double sink = 0.0;
    long double checksum = 0.0;
    auto best = std::chrono::steady_clock::duration::max();
    for (int rep = 0; rep < reps; ++rep) {
        long double chk = 0.0;
        const auto t0 = std::chrono::steady_clock::now();
        for (const std::size_t i : order) {
            const auto& fs = fss[i];
            const unsigned reg = fs.pvtRegionIndex();
            for (int ph = 0; ph < numPhases; ++ph) {
                if (!FluidSystem::phaseIsActive(ph)) {
                    continue;
                }
                const auto [b, mu] =
                    FluidSystem::inverseFormationVolumeFactorAndViscosity(fs, ph, reg);
                chk += Opm::getValue(b) + Opm::getValue(mu);
            }
        }
        const auto dt = std::chrono::steady_clock::now() - t0;
        best = std::min(best, dt);
        sink = sink + static_cast<double>(chk);
        if (rep == 0) {
            checksum = chk;
        }
    }
    (void)sink;
    return { toNsPerCell(best, order.size()), checksum };
}

// Region-batched variant C, relperm pass: SATNUM-grouped, monomorphic
// EclDefaultMaterial (no multiplexer switch in the inner loop).
template <class Evaluation>
RunResult runRelpermBatched(const Manager& manager,
                            const std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>>& fss,
                            const std::vector<std::vector<std::size_t>>& satGroups,
                            int reps)
{
    using DefaultMaterial = typename MaterialLaw::DefaultMaterial;
    using DefaultParams = typename DefaultMaterial::Params;
    using EvalArr = std::array<Evaluation, numPhases>;

    const std::size_t nc = fss.size();
    // per-cell monomorphic params, hoisted out of the timed loop
    std::vector<const DefaultParams*> params(nc, nullptr);
    for (std::size_t i = 0; i < nc; ++i) {
        params[i] = &manager.materialLawParams(i)
                        .template getRealParams<Opm::EclMultiplexerApproach::Default>();
    }

    volatile double sink = 0.0;
    long double checksum = 0.0;
    auto best = std::chrono::steady_clock::duration::max();
    std::size_t count = 0;
    for (int rep = 0; rep < reps; ++rep) {
        long double chk = 0.0;
        count = 0;
        const auto t0 = std::chrono::steady_clock::now();
        for (const auto& group : satGroups) {
            for (const std::size_t i : group) {
                const auto& fs = fss[i];
                EvalArr kr, pC;
                DefaultMaterial::relativePermeabilities(kr, *params[i], fs);
                DefaultMaterial::capillaryPressures(pC, *params[i], fs);
                for (int ph = 0; ph < numPhases; ++ph) {
                    chk += Opm::getValue(kr[ph]);
                }
                chk += Opm::getValue(pC[gasPhaseIdx]) - Opm::getValue(pC[oilPhaseIdx]);
                chk += Opm::getValue(pC[waterPhaseIdx]) - Opm::getValue(pC[oilPhaseIdx]);
                ++count;
            }
        }
        const auto dt = std::chrono::steady_clock::now() - t0;
        best = std::min(best, dt);
        sink = sink + static_cast<double>(chk);
        if (rep == 0) {
            checksum = chk;
        }
    }
    (void)sink;
    return { toNsPerCell(best, count), checksum };
}

// Variant F: the generic law stack "nested out" -- one function computing all
// six curve values, deriving each two-phase abscissa once instead of letting
// krw/krn/krg/pcnw/pcgn each re-read the fluid state and re-derive it, and
// calling the two-phase laws directly (no multiplexer dispatch). The
// representation is unchanged, so this still performs six segment searches:
// it isolates how much of the generic cost is argument re-derivation and
// dispatch rather than the searches themselves.
template <class Evaluation, class LawT = MaterialLaw, class MgrT = Manager>
RunResult runRelpermFlattened(const MgrT& manager,
                              const std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>>& fss,
                              const std::vector<std::size_t>& order,
                              int reps)
{
    using DefaultMaterial = typename LawT::DefaultMaterial;
    using DefaultParams = typename DefaultMaterial::Params;
    using GasOilLaw = typename DefaultMaterial::GasOilMaterialLaw;
    using OilWaterLaw = typename DefaultMaterial::OilWaterMaterialLaw;

    std::vector<const DefaultParams*> params(fss.size(), nullptr);
    for (std::size_t i = 0; i < fss.size(); ++i) {
        params[i] = &manager.materialLawParams(i)
                        .template getRealParams<Opm::EclMultiplexerApproach::Default>();
    }

    volatile double sink = 0.0;
    long double checksum = 0.0;
    auto best = std::chrono::steady_clock::duration::max();
    for (int rep = 0; rep < reps; ++rep) {
        long double chk = 0.0;
        const auto t0 = std::chrono::steady_clock::now();
        for (const std::size_t i : order) {
            const auto& fs = fss[i];
            const DefaultParams& p = *params[i];

            // the three independent abscissas, derived once
            const Evaluation swRaw = fs.saturation(waterPhaseIdx);
            const Evaluation sg = fs.saturation(gasPhaseIdx);
            const Scalar Swco = p.Swl();
            const Evaluation swClamped = Opm::max(Evaluation(Swco), swRaw);
            const Evaluation S = swClamped + sg;          // oil-water oil arg
            const Evaluation SoGo = 1.0 - S;              // gas-oil oil arg
            const Evaluation SgArg = 1.0 - Swco - sg;     // gas-oil gas/pc arg

            const Evaluation krw = OilWaterLaw::twoPhaseSatKrw(p.oilWaterParams(), swRaw);
            const Evaluation pcnw = OilWaterLaw::twoPhaseSatPcnw(p.oilWaterParams(), swRaw);
            const Evaluation kroOw = OilWaterLaw::twoPhaseSatKrn(p.oilWaterParams(), S);
            const Evaluation krg = GasOilLaw::twoPhaseSatKrn(p.gasOilParams(), SgArg);
            const Evaluation pcgn = GasOilLaw::twoPhaseSatPcnw(p.gasOilParams(), SgArg);
            const Evaluation kroGo = GasOilLaw::twoPhaseSatKrw(p.gasOilParams(), SoGo);

            // ECL default three-phase oil, same regularization as krn()
            Evaluation kro;
            constexpr Scalar epsilon = 1e-5;
            const Scalar Sval = Opm::getValue(S);
            if (Sval - Swco < epsilon) {
                const Evaluation kro2 = (kroOw + kroGo) / 2;
                if (Sval - Swco > epsilon / 2) {
                    const Evaluation kro1 = (sg * kroGo + (swClamped - Swco) * kroOw) / (S - Swco);
                    const Evaluation alpha = (epsilon - (S - Swco)) / (epsilon / 2);
                    kro = kro2 * alpha + kro1 * (1 - alpha);
                }
                else {
                    kro = kro2;
                }
            }
            else {
                kro = (sg * kroGo + (swClamped - Swco) * kroOw) / (S - Swco);
            }

            chk += Opm::getValue(krw) + Opm::getValue(kro) + Opm::getValue(krg);
            chk += Opm::getValue(pcgn);        // pC[gas] - pC[oil]
            chk += -Opm::getValue(pcnw);       // pC[water] - pC[oil]
        }
        const auto dt = std::chrono::steady_clock::now() - t0;
        best = std::min(best, dt);
        sink = sink + static_cast<double>(chk);
        if (rep == 0) {
            checksum = chk;
        }
    }
    (void)sink;
    return { toNsPerCell(best, order.size()), checksum };
}

// Region-batched variant C, PVT pass: PVTNUM-grouped with concrete PVT objects
// hoisted per region. Gas approach is a template parameter (WetGas or DryGas).
template <class Evaluation, Opm::GasPvtApproach gasApproach>
RunResult runPvtBatched(const std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>>& fss,
                        const std::vector<std::vector<std::size_t>>& pvtGroups,
                        int reps)
{
    // by-value copies: the concrete fused methods are (spuriously) non-const
    auto oilPvt =
        FluidSystem::oilPvt().template getRealPvt<Opm::OilPvtApproach::LiveOil>();
    auto gasPvt =
        FluidSystem::gasPvt().template getRealPvt<gasApproach>();
    auto waterPvt =
        FluidSystem::waterPvt().template getRealPvt<Opm::WaterPvtApproach::ConstantCompressibilityWater>();

    volatile double sink = 0.0;
    long double checksum = 0.0;
    auto best = std::chrono::steady_clock::duration::max();
    std::size_t count = 0;
    for (int rep = 0; rep < reps; ++rep) {
        long double chk = 0.0;
        count = 0;
        const auto t0 = std::chrono::steady_clock::now();
        for (std::size_t reg = 0; reg < pvtGroups.size(); ++reg) {
            for (const std::size_t i : pvtGroups[reg]) {
                const auto& fs = fss[i];
                const auto [bo, muo] = oilPvt.inverseFormationVolumeFactorAndViscosity(fs, reg);
                const auto [bg, mug] = gasPvt.inverseFormationVolumeFactorAndViscosity(fs, reg);
                const auto [bw, muw] = waterPvt.inverseFormationVolumeFactorAndViscosity(fs, reg);
                chk += Opm::getValue(bo) + Opm::getValue(muo)
                     + Opm::getValue(bg) + Opm::getValue(mug)
                     + Opm::getValue(bw) + Opm::getValue(muw);
                ++count;
            }
        }
        const auto dt = std::chrono::steady_clock::now() - t0;
        best = std::min(best, dt);
        sink = sink + static_cast<double>(chk);
        if (rep == 0) {
            checksum = chk;
        }
    }
    (void)sink;
    return { toNsPerCell(best, count), checksum };
}

// ---------------------------------------------------------------------------
// Variant D: baked relperm tables + analytic-derivative PVT
// ---------------------------------------------------------------------------
using Baked = Opm::Bench::BakedRelperm<MaterialLaw, Manager, Scalar>;

template <class Evaluation, class BakedT>
RunResult runRelpermBakedT(const BakedT& baked,
                           const std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>>& fss,
                           const std::vector<std::size_t>& order,
                           int reps)
{
    volatile double sink = 0.0;
    long double checksum = 0.0;
    auto best = std::chrono::steady_clock::duration::max();
    for (int rep = 0; rep < reps; ++rep) {
        long double chk = 0.0;
        const auto t0 = std::chrono::steady_clock::now();
        for (const std::size_t i : order) {
            const auto& fs = fss[i];
            std::array<Evaluation, 3> kr;
            Evaluation pcowDiff, pcgoDiff;
            baked.evaluate(i, fs.saturation(waterPhaseIdx), fs.saturation(gasPhaseIdx),
                           kr, pcowDiff, pcgoDiff);
            for (int ph = 0; ph < numPhases; ++ph) {
                chk += Opm::getValue(kr[ph]);
            }
            chk += Opm::getValue(pcgoDiff) + Opm::getValue(pcowDiff);
        }
        const auto dt = std::chrono::steady_clock::now() - t0;
        best = std::min(best, dt);
        sink = sink + static_cast<double>(chk);
        if (rep == 0) { checksum = chk; }
    }
    (void)sink;
    return { toNsPerCell(best, order.size()), checksum };
}

// Replay a saturation history so both two-phase systems record a turning
// point: one pass at higher Sg (gas-oil MDC) and one at higher So (oil-water
// MDC). update() keeps the minimum of its argument, so the two passes are
// independent — afterwards the original states lie on the imbibition side of
// both curves.
template <class MgrT>
void driveIntoImbibition(MgrT& manager, const std::vector<ResolvedCell>& cells, std::size_t nc)
{
    for (std::size_t c = 0; c < nc; ++c) {
        auto hi = makeFluidState<Scalar>(cells[c], 0);
        const Scalar sw = cells[c].sw, sg = cells[c].sg;

        const Scalar sgHi = std::min(Scalar(1.0 - sw - 1e-3), Scalar(sg + 0.2));
        hi.setSaturation(waterPhaseIdx, sw);
        hi.setSaturation(gasPhaseIdx, sgHi);
        hi.setSaturation(oilPhaseIdx, 1.0 - sw - sgHi);
        manager.updateHysteresis(hi, static_cast<unsigned>(c));

        const Scalar swLo = std::max(Scalar(0.0), Scalar(sw - 0.1));
        const Scalar sgLo = std::max(Scalar(0.0), Scalar(sg - 0.1));
        hi.setSaturation(waterPhaseIdx, swLo);
        hi.setSaturation(gasPhaseIdx, sgLo);
        hi.setSaturation(oilPhaseIdx, 1.0 - swLo - sgLo);
        manager.updateHysteresis(hi, static_cast<unsigned>(c));
    }
}

// Cell-by-cell agreement of the baked path against the true law, values and
// derivatives, plus how many cells actually sit on an imbibition branch.
template <class Evaluation, class LawT, class MgrT, class BakedT>
void compareBakedVsGeneric(const MgrT& manager,
                           const BakedT& baked,
                           const std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>>& fss,
                           std::size_t nc)
{
    using EvalArr = std::array<Evaluation, numPhases>;
    double maxKr = 0.0, maxPc = 0.0, maxDeriv = 0.0;
    std::size_t revGo = 0, revOw = 0;
    std::size_t worstCell = 0, nDerivBad = 0;
    double worstDen = 1e30;
    int worstQty = -1, worstK = -1;

    for (std::size_t c = 0; c < nc; ++c) {
        const auto& fs = fss[c];
        const auto& mp = manager.materialLawParams(c);
        const auto& dp = mp.template getRealParams<Opm::EclMultiplexerApproach::Default>();

        const double sw = Opm::getValue(fs.saturation(waterPhaseIdx));
        const double sg = Opm::getValue(fs.saturation(gasPhaseIdx));
        if (1.0 - dp.Swl() - sg > dp.gasOilParams().krnSwMdc())  { ++revGo; }
        if (sw + sg > dp.oilWaterParams().krnSwMdc())            { ++revOw; }

        EvalArr kr, pC;
        LawT::relativePermeabilities(kr, mp, fs);
        LawT::capillaryPressures(pC, mp, fs);

        std::array<Evaluation, 3> krB;
        Evaluation pcowB, pcgoB;
        baked.evaluate(c, fs.saturation(waterPhaseIdx), fs.saturation(gasPhaseIdx),
                       krB, pcowB, pcgoB);

        const Evaluation pcow = pC[waterPhaseIdx] - pC[oilPhaseIdx];
        const Evaluation pcgo = pC[gasPhaseIdx] - pC[oilPhaseIdx];
        auto track = [](double& acc, const Evaluation& a, const Evaluation& b) {
            acc = std::max(acc, std::abs(Opm::getValue(a) - Opm::getValue(b)));
        };
        auto trackD = [&](const Evaluation& a, const Evaluation& b, int qty) {
            for (int k = 0; k < static_cast<int>(a.size()) - 1; ++k) {
                const double d = std::abs(a.derivative(k) - b.derivative(k));
                if (d > 1e-3) {
                    ++nDerivBad;
                    // the three-phase kro combine divides by (S - Swco);
                    // record how near-singular the offenders are
                    worstDen = std::min(worstDen, std::max(sw, dp.Swl()) + sg - dp.Swl());
                }
                if (d > maxDeriv) { maxDeriv = d; worstCell = c; worstQty = qty; worstK = k; }
            }
        };
        for (int ph = 0; ph < numPhases; ++ph) {
            track(maxKr, kr[ph], krB[ph]);
            trackD(kr[ph], krB[ph], ph);
        }
        track(maxPc, pcow, pcowB);
        track(maxPc, pcgo, pcgoB);
        trackD(pcow, pcowB, 3);
        trackD(pcgo, pcgoB, 4);
    }

    std::printf("imbibition check: %zu/%zu cells reversed gas-oil, %zu/%zu oil-water; "
                "max |dkr| %.2e  max |dpc| %.2e Pa  max |dderiv| %.2e "
                "(%zu deriv comps >1e-3, min |S-Swco| %.2e; worst cell %zu qty %d d/dx%d)\n",
                revGo, nc, revOw, nc, maxKr, maxPc, maxDeriv,
                nDerivBad, nDerivBad ? worstDen : 0.0, worstCell, worstQty, worstK);
    if (worstQty >= 0) {
        const auto& mp = manager.materialLawParams(worstCell);
        const auto& dp = mp.template getRealParams<Opm::EclMultiplexerApproach::Default>();
        std::printf("  worst cell: sw %.6f sg %.6f Swl %.6f | go krnSwMdc %.6f delta %.3e"
                    " | ow krnSwMdc %.6f delta %.3e\n",
                    Opm::getValue(fss[worstCell].saturation(waterPhaseIdx)),
                    Opm::getValue(fss[worstCell].saturation(gasPhaseIdx)),
                    dp.Swl(),
                    dp.gasOilParams().krnSwMdc(), dp.gasOilParams().deltaSwImbKrn(),
                    dp.oilWaterParams().krnSwMdc(), dp.oilWaterParams().deltaSwImbKrn());
    }
}

template <class Evaluation, bool uniform = false>
RunResult runRelpermBaked(const Baked& baked,
                          const std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>>& fss,
                          const std::vector<std::size_t>& order,
                          int reps)
{
    volatile double sink = 0.0;
    long double checksum = 0.0;
    auto best = std::chrono::steady_clock::duration::max();
    for (int rep = 0; rep < reps; ++rep) {
        long double chk = 0.0;
        const auto t0 = std::chrono::steady_clock::now();
        for (const std::size_t i : order) {
            const auto& fs = fss[i];
            std::array<Evaluation, 3> kr;
            Evaluation pcowDiff, pcgoDiff;
            if constexpr (uniform) {
                baked.evaluateUniform(i, fs.saturation(waterPhaseIdx), fs.saturation(gasPhaseIdx),
                                      kr, pcowDiff, pcgoDiff);
            }
            else {
                baked.evaluate(i, fs.saturation(waterPhaseIdx), fs.saturation(gasPhaseIdx),
                               kr, pcowDiff, pcgoDiff);
            }
            for (int ph = 0; ph < numPhases; ++ph) {
                chk += Opm::getValue(kr[ph]);
            }
            chk += Opm::getValue(pcgoDiff) + Opm::getValue(pcowDiff);
        }
        const auto dt = std::chrono::steady_clock::now() - t0;
        best = std::min(best, dt);
        sink = sink + static_cast<double>(chk);
        if (rep == 0) {
            checksum = chk;
        }
    }
    (void)sink;
    return { toNsPerCell(best, order.size()), checksum };
}

// 1D piecewise-linear eval with analytic derivative construction.
template <class Evaluation>
Evaluation eval1dAnalytic(const Opm::Tabulated1DFunction<Scalar>& tab, const Evaluation& xe)
{
    const auto& xs = tab.xValues();
    const auto& ys = tab.yValues();
    const double xd = Opm::getValue(xe);
    const std::size_t n = xs.size();
    std::size_t i = std::upper_bound(xs.begin(), xs.end(), xd) - xs.begin();
    i = std::min(std::max<std::size_t>(i, 1), n - 1) - 1;
    const double slope = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]);
    if constexpr (std::is_same_v<Evaluation, Scalar>) {
        return ys[i] + slope * (xd - xs[i]);
    }
    else {
        Evaluation r;
        r.setValue(ys[i] + slope * (xd - xs[i]));
        for (int k = 0; k < xe.size(); ++k) {
            r.setDerivative(k, slope * xe.derivative(k));
        }
        return r;
    }
}

// 2D eval on UniformXTabulated2DFunction with analytic partials.
template <class Evaluation>
Evaluation eval2dAnalytic(const Opm::UniformXTabulated2DFunction<Scalar>& tab,
                          const Evaluation& xe, const Evaluation& ye)
{
    unsigned ii, j1, j2;
    double alpha, beta1, beta2;
    tab.findPoints(ii, j1, j2, alpha, beta1, beta2,
                   Opm::getValue(xe), Opm::getValue(ye), /*extrapolate=*/true);
    const double a0 = tab.valueAt(ii, j1);
    const double a1 = tab.valueAt(ii, j1 + 1);
    const double b0 = tab.valueAt(ii + 1, j2);
    const double b1 = tab.valueAt(ii + 1, j2 + 1);
    const double s1v = a0 * (1.0 - beta1) + a1 * beta1;
    const double s2v = b0 * (1.0 - beta2) + b1 * beta2;
    const double v = s1v * (1.0 - alpha) + s2v * alpha;
    if constexpr (std::is_same_v<Evaluation, Scalar>) {
        return v;
    }
    else {
        const double dvdx = (s2v - s1v) / (tab.iToX(ii + 1) - tab.iToX(ii));
        const double dvdy =
            (1.0 - alpha) * (a1 - a0) / (tab.yAt(ii, j1 + 1) - tab.yAt(ii, j1))
            + alpha * (b1 - b0) / (tab.yAt(ii + 1, j2 + 1) - tab.yAt(ii + 1, j2));
        Evaluation r;
        r.setValue(v);
        for (int k = 0; k < xe.size(); ++k) {
            r.setDerivative(k, dvdx * xe.derivative(k) + dvdy * ye.derivative(k));
        }
        return r;
    }
}

// Analytic-derivative PVT (LiveOil + Wet/DryGas tables, CCWater closed form),
// mirroring the branch logic of the fused inverseFormationVolumeFactorAndViscosity.
template <class Evaluation, Opm::GasPvtApproach gasApproach>
RunResult runPvtAnalytic(const std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>>& fss,
                         const std::vector<std::vector<std::size_t>>& pvtGroups,
                         int reps)
{
    constexpr bool dryGas = (gasApproach == Opm::GasPvtApproach::DryGas);
    const auto& oilPvt =
        FluidSystem::oilPvt().template getRealPvt<Opm::OilPvtApproach::LiveOil>();
    const auto& gasPvt =
        FluidSystem::gasPvt().template getRealPvt<gasApproach>();
    auto waterPvt =
        FluidSystem::waterPvt().template getRealPvt<Opm::WaterPvtApproach::ConstantCompressibilityWater>();

    volatile double sink = 0.0;
    long double checksum = 0.0;
    auto best = std::chrono::steady_clock::duration::max();
    std::size_t count = 0;
    for (int rep = 0; rep < reps; ++rep) {
        long double chk = 0.0;
        count = 0;
        const auto t0 = std::chrono::steady_clock::now();
        for (std::size_t reg = 0; reg < pvtGroups.size(); ++reg) {
            const auto& rsSatTab = oilPvt.saturatedGasDissolutionFactorTable()[reg];
            const auto& oilBSat = oilPvt.inverseSaturatedOilBTable()[reg];
            const auto& oilBMuSat = oilPvt.inverseSaturatedOilBMuTable()[reg];
            const auto& oilB = oilPvt.inverseOilBTable()[reg];
            const auto& oilBMu = oilPvt.inverseOilBMuTable()[reg];
            const auto& gasB = gasPvt.inverseGasB()[reg];
            const auto& gasBMu = gasPvt.inverseGasBMu()[reg];

            for (const std::size_t i : pvtGroups[reg]) {
                const auto& fs = fss[i];
                const Evaluation& po = fs.pressure(oilPhaseIdx);
                const Evaluation& pg = fs.pressure(gasPhaseIdx);
                const Evaluation& pw = fs.pressure(waterPhaseIdx);
                const Evaluation& rs = fs.Rs();
                const Evaluation& rv = fs.Rv();

                // oil
                {
                    const Evaluation rsSat = eval1dAnalytic(rsSatTab, po);
                    const bool sat = (Opm::getValue(fs.saturation(gasPhaseIdx)) > 0.0)
                        && (Opm::getValue(rs) >= (1.0 - 1e-10) * Opm::getValue(rsSat));
                    Evaluation b, ibmu;
                    if (sat) {
                        b = eval1dAnalytic(oilBSat, po);
                        ibmu = eval1dAnalytic(oilBMuSat, po);
                    }
                    else {
                        b = eval2dAnalytic(oilB, rs, po);
                        ibmu = eval2dAnalytic(oilBMu, rs, po);
                    }
                    const Evaluation mu = b / ibmu;
                    chk += Opm::getValue(b) + Opm::getValue(mu);
                }
                // gas
                if constexpr (dryGas) {
                    const Evaluation b = eval1dAnalytic(gasB, pg);
                    const Evaluation ibmu = eval1dAnalytic(gasBMu, pg);
                    const Evaluation mu = b / ibmu;
                    chk += Opm::getValue(b) + Opm::getValue(mu);
                }
                else {
                    const auto& rvSatTab = gasPvt.saturatedOilVaporizationFactorTable()[reg];
                    const auto& gasBSat = gasPvt.inverseSaturatedGasB()[reg];
                    const auto& gasBMuSat = gasPvt.inverseSaturatedGasBMu()[reg];
                    const Evaluation rvSat = eval1dAnalytic(rvSatTab, pg);
                    const bool sat = (Opm::getValue(fs.saturation(oilPhaseIdx)) > 0.0)
                        && (Opm::getValue(rv) >= (1.0 - 1e-10) * Opm::getValue(rvSat));
                    Evaluation b, ibmu;
                    if (sat) {
                        b = eval1dAnalytic(gasBSat, pg);
                        ibmu = eval1dAnalytic(gasBMuSat, pg);
                    }
                    else {
                        b = eval2dAnalytic(gasB, pg, rv);
                        ibmu = eval2dAnalytic(gasBMu, pg, rv);
                    }
                    const Evaluation mu = b / ibmu;
                    chk += Opm::getValue(b) + Opm::getValue(mu);
                }
                // water: closed-form, no tables
                {
                    Evaluation bw, muw;
                    waterPvt.inverseBAndMu(bw, muw, static_cast<unsigned>(reg), pw);
                    chk += Opm::getValue(bw) + Opm::getValue(muw);
                }
                ++count;
            }
        }
        const auto dt = std::chrono::steady_clock::now() - t0;
        best = std::min(best, dt);
        sink = sink + static_cast<double>(chk);
        if (rep == 0) {
            checksum = chk;
        }
    }
    (void)sink;
    return { toNsPerCell(best, count), checksum };
}

// The IQ's actual PVT sequence: first the composition step evaluates RsSat /
// RvSat to set Rs/Rv, then the b+mu step needs the same saturated value to
// pick its branch. Master recomputes the segment search and the table eval;
// `shareIndex` carries them over instead (plan §2a). Everything else is
// identical, so the delta is exactly the redundant search+eval.
template <class Evaluation, Opm::GasPvtApproach gasApproach, bool shareIndex,
          bool shortCircuit = false>
RunResult runPvtIQLike(const std::vector<Opm::BlackOilFluidState<Evaluation, FluidSystem>>& fss,
                       const std::vector<std::vector<std::size_t>>& pvtGroups,
                       int reps)
{
    constexpr bool dryGas = (gasApproach == Opm::GasPvtApproach::DryGas);
    const auto& oilPvt =
        FluidSystem::oilPvt().template getRealPvt<Opm::OilPvtApproach::LiveOil>();
    const auto& gasPvt =
        FluidSystem::gasPvt().template getRealPvt<gasApproach>();
    auto waterPvt =
        FluidSystem::waterPvt().template getRealPvt<Opm::WaterPvtApproach::ConstantCompressibilityWater>();

    volatile double sink = 0.0;
    long double checksum = 0.0;
    auto best = std::chrono::steady_clock::duration::max();
    std::size_t count = 0;
    for (int rep = 0; rep < reps; ++rep) {
        long double chk = 0.0;
        count = 0;
        const auto t0 = std::chrono::steady_clock::now();
        for (std::size_t reg = 0; reg < pvtGroups.size(); ++reg) {
            const auto& rsSatTab = oilPvt.saturatedGasDissolutionFactorTable()[reg];
            const auto& oilBSat = oilPvt.inverseSaturatedOilBTable()[reg];
            const auto& oilBMuSat = oilPvt.inverseSaturatedOilBMuTable()[reg];
            const auto& oilB = oilPvt.inverseOilBTable()[reg];
            const auto& oilBMu = oilPvt.inverseOilBMuTable()[reg];
            const auto& gasB = gasPvt.inverseGasB()[reg];
            const auto& gasBMu = gasPvt.inverseGasBMu()[reg];

            for (const std::size_t i : pvtGroups[reg]) {
                const auto& fs = fss[i];
                const Evaluation& po = fs.pressure(oilPhaseIdx);
                const Evaluation& pg = fs.pressure(gasPhaseIdx);
                const Evaluation& pw = fs.pressure(waterPhaseIdx);
                const Evaluation& rs = fs.Rs();
                const Evaluation& rv = fs.Rv();

                // --- oil: composition step then b+mu ---
                // Mirrors the IQ: with GasMeaning::Rs (no free gas) the
                // composition step reads Rs from the primary variables and
                // never evaluates RsSat, yet the b/mu step evaluates it only
                // to have the flag discard it (saturation(gas) > 0 fails).
                {
                    const bool hasGas = Opm::getValue(fs.saturation(gasPhaseIdx)) > 0.0;
                    Evaluation rsSatComp;
                    Opm::SegmentIndex segIdx{0};
                    bool haveComp = false;
                    if (hasGas) {
                        segIdx = rsSatTab.findSegmentIndex(Opm::getValue(po), true);
                        rsSatComp = rsSatTab.eval(po, segIdx);
                        chk += Opm::getValue(rsSatComp);       // stands for setRs()
                        haveComp = true;
                    }

                    bool sat = false;
                    const bool needSaturatedFactor = hasGas || !shortCircuit;
                    if (needSaturatedFactor) {
                        Evaluation rsSatUse;
                        if (shareIndex && haveComp) {
                            rsSatUse = rsSatComp;              // carried over
                        }
                        else {
                            segIdx = rsSatTab.findSegmentIndex(Opm::getValue(po), true);
                            rsSatUse = rsSatTab.eval(po, segIdx);
                        }
                        sat = hasGas
                            && (Opm::getValue(rs) >= (1.0 - 1e-10) * Opm::getValue(rsSatUse));
                    }

                    Evaluation b, ibmu;
                    if (sat) {
                        b = oilBSat.eval(po, segIdx);
                        ibmu = oilBMuSat.eval(po, segIdx);
                    }
                    else {
                        unsigned ii, j1, j2;
                        Evaluation alpha, beta1, beta2;
                        oilB.findPoints(ii, j1, j2, alpha, beta1, beta2, rs, po, true);
                        b = oilB.eval(ii, j1, j2, alpha, beta1, beta2);
                        ibmu = oilBMu.eval(ii, j1, j2, alpha, beta1, beta2);
                    }
                    chk += Opm::getValue(b) + Opm::getValue(b / ibmu);
                }
                // --- gas ---
                if constexpr (dryGas) {
                    const Evaluation b = gasB.eval(pg, true);
                    const Evaluation ibmu = gasBMu.eval(pg, true);
                    chk += Opm::getValue(b) + Opm::getValue(b / ibmu);
                }
                else {
                    const auto& rvSatTab = gasPvt.saturatedOilVaporizationFactorTable()[reg];
                    const auto& gasBSat = gasPvt.inverseSaturatedGasB()[reg];
                    const auto& gasBMuSat = gasPvt.inverseSaturatedGasBMu()[reg];
                    auto segIdx = rvSatTab.findSegmentIndex(Opm::getValue(pg), true);
                    const Evaluation rvSatComp = rvSatTab.eval(pg, segIdx);
                    chk += Opm::getValue(rvSatComp);          // stands for setRv()

                    Evaluation rvSatUse;
                    if constexpr (shareIndex) {
                        rvSatUse = rvSatComp;
                    }
                    else {
                        segIdx = rvSatTab.findSegmentIndex(Opm::getValue(pg), true);
                        rvSatUse = rvSatTab.eval(pg, segIdx);
                    }
                    const bool sat = (Opm::getValue(fs.saturation(oilPhaseIdx)) > 0.0)
                        && (Opm::getValue(rv) >= (1.0 - 1e-10) * Opm::getValue(rvSatUse));
                    Evaluation b, ibmu;
                    if (sat) {
                        b = gasBSat.eval(pg, segIdx);
                        ibmu = gasBMuSat.eval(pg, segIdx);
                    }
                    else {
                        unsigned ii, j1, j2;
                        Evaluation alpha, beta1, beta2;
                        gasB.findPoints(ii, j1, j2, alpha, beta1, beta2, pg, rv, true);
                        b = gasB.eval(ii, j1, j2, alpha, beta1, beta2);
                        ibmu = gasBMu.eval(ii, j1, j2, alpha, beta1, beta2);
                    }
                    chk += Opm::getValue(b) + Opm::getValue(b / ibmu);
                }
                // --- water ---
                {
                    Evaluation bw, muw;
                    waterPvt.inverseBAndMu(bw, muw, static_cast<unsigned>(reg), pw);
                    chk += Opm::getValue(bw) + Opm::getValue(muw);
                }
                ++count;
            }
        }
        const auto dt = std::chrono::steady_clock::now() - t0;
        best = std::min(best, dt);
        sink = sink + static_cast<double>(chk);
        if (rep == 0) {
            checksum = chk;
        }
    }
    (void)sink;
    return { toNsPerCell(best, count), checksum };
}

void printRow(const std::string& what, const std::string& detail, const RunResult& r)
{
    std::printf("%-34s %-10s %10.1f ns/cell   chk %.10Le\n",
                what.c_str(), detail.c_str(), r.nsPerCell, r.checksum);
}

} // namespace

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    const Options opts = parseArgs(argc, argv);

    // --- deck / state ------------------------------------------------------
    Opm::ParseContext parseContext({{Opm::ParseContext::PARSE_RANDOM_SLASH,
                                     Opm::InputErrorAction::IGNORE}});
    Opm::ErrorGuard errors;
    Opm::Parser parser;
    const auto deck = parser.parseFile(opts.deck, parseContext, errors);
    const auto python = std::make_shared<Opm::Python>();
    const Opm::EclipseState eclState(deck);
    const Opm::Schedule schedule(deck, eclState, python);

    FluidSystem::initFromState(eclState, schedule);

    std::size_t nc = eclState.getInputGrid().getNumActive();
    if (opts.n > 0 && opts.n < nc) {
        nc = opts.n;
    }

    Manager manager;
    manager.initFromState(eclState);
    {
        // LGR-aware init signature; identity lookups (pattern from the
        // existing material tests).
        std::function<std::vector<int>(const Opm::FieldPropsManager&, const std::string&, bool)>
            lookup = [](const Opm::FieldPropsManager& fp, const std::string& prop, bool translate) {
                std::vector<int> dest;
                const auto& raw = fp.get_int(prop);
                dest.resize(raw.size());
                for (std::size_t i = 0; i < raw.size(); ++i) {
                    dest[i] = raw[i] - static_cast<int>(translate);
                }
                return dest;
            };
        std::function<unsigned(unsigned)> identity = [](unsigned i) { return i; };
        manager.initParamsForElements(eclState, eclState.getInputGrid().getNumActive(),
                                      lookup, identity);
    }
    if (opts.checkOnly) {
        // Authoritative coverage report: same parser, same approach selection
        // and same guard a simulator would hit. Uses hysteresis-enabled traits
        // because that is the configuration a deck with EHYSTR would run.
        using HystTraits = Opm::ThreePhaseMaterialTraits<Scalar,
                                                         FluidSystem::waterPhaseIdx,
                                                         FluidSystem::oilPhaseIdx,
                                                         FluidSystem::gasPhaseIdx,
                                                         true, true>;
        using HystManager = Opm::EclMaterialLaw::Manager<HystTraits>;
        HystManager hm;
        hm.initFromState(eclState);
        std::function<std::vector<int>(const Opm::FieldPropsManager&, const std::string&, bool)>
            lk = [](const Opm::FieldPropsManager& fp, const std::string& prop, bool translate) {
                const auto& raw = fp.get_int(prop);
                std::vector<int> dest(raw.size());
                for (std::size_t i = 0; i < raw.size(); ++i) {
                    dest[i] = raw[i] - static_cast<int>(translate);
                }
                return dest;
            };
        std::function<unsigned(unsigned)> id = [](unsigned i) { return i; };
        hm.initParamsForElements(eclState, eclState.getInputGrid().getNumActive(), lk, id);
        using HystBaked = Opm::Bench::BakedRelperm<typename HystManager::MaterialLaw,
                                                   HystManager, Scalar>;
        const auto why = HystBaked::unsupportedReason(hm);
        static const char* approachName[] = { "Default", "Stone1", "Stone2",
                                              "TwoPhase", "OnePhase" };
        std::printf("%-40s %-9s hyst=%d krModel=%2d pcModel=%2d wag=%d  %s\n",
                    opts.deck.c_str(),
                    approachName[static_cast<int>(hm.threePhaseApproach())],
                    int(hm.enableHysteresis()),
                    hm.hysteresisConfig().krHysteresisModel(),
                    hm.hysteresisConfig().pcHysteresisModel(),
                    int(hm.hysteresisConfig().enableWagHysteresis()),
                    why.empty() ? "SUPPORTED" : ("unsupported: " + why).c_str());
        return 0;
    }
    if (manager.hasDirectionalRelperms() || manager.hasDirectionalImbnum()) {
        std::cerr << "directional relperms not supported by this benchmark\n";
        return 1;
    }

    // region arrays (0-based); never sorted in place — permutations only
    std::vector<int> pvtnum(nc, 0);
    if (eclState.fieldProps().has_int("PVTNUM")) {
        const auto& raw = eclState.fieldProps().get_int("PVTNUM");
        for (std::size_t i = 0; i < nc; ++i) {
            pvtnum[i] = raw[i] - 1;
        }
    }
    std::vector<int> satnum(nc, 0);
    for (std::size_t i = 0; i < nc; ++i) {
        satnum[i] = static_cast<int>(manager.satnumRegionIdx(i));
    }

    std::vector<Scalar> poro(nc, 0.25);
    if (eclState.fieldProps().has_double("PORO")) {
        const auto& raw = eclState.fieldProps().get_double("PORO");
        for (std::size_t i = 0; i < nc; ++i) {
            poro[i] = raw[i];
        }
    }

    // orders: natural, region-sorted, shuffled
    std::vector<std::size_t> natural(nc);
    std::iota(natural.begin(), natural.end(), 0);
    std::vector<std::size_t> sorted = natural;
    std::stable_sort(sorted.begin(), sorted.end(),
                     [&](std::size_t a, std::size_t b) {
                         return std::make_pair(pvtnum[a], satnum[a])
                              < std::make_pair(pvtnum[b], satnum[b]);
                     });
    std::vector<std::size_t> shuffled = natural;
    std::shuffle(shuffled.begin(), shuffled.end(), std::mt19937(opts.seed + 1));

    // region groups (CSR-style buckets) for variant C
    const auto makeGroups = [nc](const std::vector<int>& reg) {
        const int nreg = *std::max_element(reg.begin(), reg.begin() + nc) + 1;
        std::vector<std::vector<std::size_t>> groups(nreg);
        for (std::size_t i = 0; i < nc; ++i) {
            groups[reg[i]].push_back(i);
        }
        return groups;
    };
    const auto satGroups = makeGroups(satnum);
    const auto pvtGroups = makeGroups(pvtnum);

    const auto states = fabricateStates(nc, pvtnum, opts.undersatFrac, opts.seed);

    // rock constants (documented stand-ins; identical for every variant)
    const Scalar rockComp = 4.0e-10;   // 1/Pa
    const Scalar rockRefP = 280.0e5;   // Pa

    std::printf("deck %s  nc %zu  undersat-frac %.2f  seed %u  reps %d\n",
                opts.deck.c_str(), nc, opts.undersatFrac, opts.seed, opts.reps);
    std::printf("pvt approaches: oil=%d gas=%d water=%d  (LiveOil=%d WetGas=%d CCWater=%d)\n",
                static_cast<int>(FluidSystem::oilPvt().approach()),
                static_cast<int>(FluidSystem::gasPvt().approach()),
                static_cast<int>(FluidSystem::waterPvt().approach()),
                static_cast<int>(Opm::OilPvtApproach::LiveOil),
                static_cast<int>(Opm::GasPvtApproach::WetGas),
                static_cast<int>(Opm::WaterPvtApproach::ConstantCompressibilityWater));

    const bool oilWaterOk =
        FluidSystem::oilPvt().approach() == Opm::OilPvtApproach::LiveOil
        && FluidSystem::waterPvt().approach() == Opm::WaterPvtApproach::ConstantCompressibilityWater;
    const bool wetGas = FluidSystem::gasPvt().approach() == Opm::GasPvtApproach::WetGas;
    const bool dryGas = FluidSystem::gasPvt().approach() == Opm::GasPvtApproach::DryGas;
    const bool nornePvtCombo = oilWaterOk && wetGas;       // fast IQ hardcodes WetGas
    const bool concretePvtOk = oilWaterOk && (wetGas || dryGas);
    const auto approach0 = manager.materialLawParams(0).approach();
    const bool defaultMultiplexer = approach0 == Opm::EclMultiplexerApproach::Default;
    // the baked path now also covers the Stone models (two axes, three columns)
    const bool bakeableApproach = defaultMultiplexer
        || approach0 == Opm::EclMultiplexerApproach::Stone1
        || approach0 == Opm::EclMultiplexerApproach::Stone2;

    // --- variant A: master IQ ---------------------------------------------
    Opm::Bench::StubIQProblem<MasterTag> problemA(manager, poro, rockComp, rockRefP,
                                                  benchTemperature);
    const auto pvsA = makePrimaryVariables<MasterTag>(states, pvtnum);
    std::printf("\n== full IQ update (TypeTag Evaluation) ==\n");
    printRow("A master IQ", "natural", runIQVariant<MasterTag>(problemA, pvsA, natural, opts.reps));
    printRow("A master IQ", "sorted", runIQVariant<MasterTag>(problemA, pvsA, sorted, opts.reps));
    printRow("A master IQ", "shuffled", runIQVariant<MasterTag>(problemA, pvsA, shuffled, opts.reps));

    // --- variant B: fast IQ (Norne PVT combo only) -------------------------
    if (nornePvtCombo) {
        Opm::Bench::StubIQProblem<FastTag> problemB(manager, poro, rockComp, rockRefP,
                                                    benchTemperature);
        const auto pvsB = makePrimaryVariables<FastTag>(states, pvtnum);
        printRow("B fast IQ", "natural", runIQVariant<FastTag>(problemB, pvsB, natural, opts.reps));
        printRow("B fast IQ", "sorted", runIQVariant<FastTag>(problemB, pvsB, sorted, opts.reps));
        printRow("B fast IQ", "shuffled", runIQVariant<FastTag>(problemB, pvsB, shuffled, opts.reps));
    }
    else {
        std::printf("B fast IQ: skipped (needs LiveOil+WetGas+CCWater)\n");
    }

    // --- resolved states for the decomposed loops --------------------------
    const auto cells = resolveStates<MasterTag>(problemA, pvsA);

    // --- sub-part split (natural order) ------------------------------------
    using Eval3 = Opm::DenseAd::Evaluation<Scalar, 3>;
    using Eval6 = Opm::DenseAd::Evaluation<Scalar, 6>;
    const auto fssD = buildFluidStates<Scalar>(cells, pvtnum);
    const auto fss3 = buildFluidStates<Eval3>(cells, pvtnum);
    const auto fss6 = buildFluidStates<Eval6>(cells, pvtnum);
    std::printf("\n== sub-parts, generic per-cell calls (natural order) ==\n");
    printRow("relperm+pc", "double", runRelperm<Scalar>(manager, fssD, natural, opts.reps));
    printRow("relperm+pc", "eval3", runRelperm<Eval3>(manager, fss3, natural, opts.reps));
    printRow("relperm+pc", "eval6", runRelperm<Eval6>(manager, fss6, natural, opts.reps));
    printRow("relperm+pc pinned params", "eval3",
             runRelperm<Eval3>(manager, fss3, natural, opts.reps, true));
    // Footprint sweep: identical hysteresis-free computation, per-cell params
    // padded out. Reads off what shrinking the 1440 B hysteresis params would
    // recover, without having to refactor opm-common first.
    {   // Is master's per-cell allocation actually scattered, or does the
        // allocator hand out a tight run because they are made back-to-back?
        std::vector<std::uintptr_t> a(nc);
        for (std::size_t i = 0; i < nc; ++i) {
            a[i] = reinterpret_cast<std::uintptr_t>(
                &manager.materialLawParams(i)
                     .template getRealParams<Opm::EclMultiplexerApproach::Default>());
        }
        std::size_t adjacent = 0;
        double meanStride = 0.0;
        for (std::size_t i = 1; i < nc; ++i) {
            const auto d = a[i] > a[i - 1] ? a[i] - a[i - 1] : a[i - 1] - a[i];
            meanStride += double(d);
            if (d <= 512) { ++adjacent; }
        }
        auto [lo, hi] = std::minmax_element(a.begin(), a.end());
        std::printf("  params layout: %.0f%% of consecutive cells within 512 B, "
                    "mean |stride| %.0f B, span %.1f MB for %zu x %zu B\n",
                    100.0 * double(adjacent) / double(nc - 1), meanStride / double(nc - 1),
                    double(*hi - *lo) / (1024.0 * 1024.0), nc,
                    sizeof(typename Manager::MaterialLaw::Params));
    }
    printRow("contiguous, native 304 B/cell", "eval3",
             runRelpermPadded<Eval3, 304>(manager, fss3, natural, opts.reps));
    printRow("contiguous, padded to  620 B", "eval3",
             runRelpermPadded<Eval3, 620>(manager, fss3, natural, opts.reps));
    printRow("contiguous, padded to  760 B", "eval3",
             runRelpermPadded<Eval3, 760>(manager, fss3, natural, opts.reps));
    printRow("contiguous, padded to 1440 B", "eval3",
             runRelpermPadded<Eval3, 1440>(manager, fss3, natural, opts.reps));
    printRow("pvt b+mu", "double", runPvt<Scalar>(fssD, natural, opts.reps));
    printRow("pvt b+mu", "eval3", runPvt<Eval3>(fss3, natural, opts.reps));
    printRow("pvt b+mu", "eval6", runPvt<Eval6>(fss6, natural, opts.reps));

    // --- order sensitivity (eval3) -----------------------------------------
    std::printf("\n== order sensitivity of generic sub-parts (eval3) ==\n");
    printRow("relperm+pc", "sorted", runRelperm<Eval3>(manager, fss3, sorted, opts.reps));
    printRow("relperm+pc", "shuffled", runRelperm<Eval3>(manager, fss3, shuffled, opts.reps));
    printRow("pvt b+mu", "sorted", runPvt<Eval3>(fss3, sorted, opts.reps));
    printRow("pvt b+mu", "shuffled", runPvt<Eval3>(fss3, shuffled, opts.reps));

    // --- variant C: region-batched, monomorphic inner loops ----------------
    std::printf("\n== region-batched passes (eval3) ==\n");
    if (defaultMultiplexer) {
        printRow("C relperm SATNUM-batched", "eval3",
                 runRelpermBatched<Eval3>(manager, fss3, satGroups, opts.reps));
        printRow("F relperm flattened (6 searches)", "double",
                 runRelpermFlattened<Scalar>(manager, fssD, natural, opts.reps));
        printRow("F relperm flattened (6 searches)", "eval3",
                 runRelpermFlattened<Eval3>(manager, fss3, natural, opts.reps));
    }
    else {
        std::printf("C relperm: skipped (multiplexer approach not Default)\n");
    }
    if (concretePvtOk) {
        constexpr auto WG = Opm::GasPvtApproach::WetGas;
        constexpr auto DG = Opm::GasPvtApproach::DryGas;
        printRow("C pvt PVTNUM-batched", "eval3",
                 wetGas ? runPvtBatched<Eval3, WG>(fss3, pvtGroups, opts.reps)
                        : runPvtBatched<Eval3, DG>(fss3, pvtGroups, opts.reps));
    }
    else {
        std::printf("C pvt: skipped (needs LiveOil+Wet/DryGas+CCWater)\n");
    }

    // --- variant D: baked relperm + analytic-derivative PVT ----------------
    std::printf("\n== variant D: baked tables + analytic derivatives ==\n");
    if (defaultMultiplexer) {
        const auto tb0 = std::chrono::steady_clock::now();
        Baked baked;
        // force: this manager runs hysteresis-free traits on purpose, so the
        // rows below stay comparable with the pre-hysteresis measurements.
        if (const auto why = Baked::unsupportedReason(manager); !why.empty()) {
            std::printf("note: baked guard would reject this build (%s); forced for comparison\n",
                        why.c_str());
        }
        baked.build(manager, nc, /*force=*/true);
        const double bakeSec =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - tb0).count();
        std::printf("bake: %zu unique tables, %zu nodes total (%.1f KB), "
                    "max validation err %.2e, build %.2f s\n",
                    baked.numTables(), baked.numNodes(),
                    baked.numNodes() * 5 * 8 / 1024.0, baked.maxBakeError(), bakeSec);
        printRow("D relperm baked", "double", runRelpermBaked<Scalar>(baked, fssD, natural, opts.reps));
        printRow("D relperm baked", "eval3", runRelpermBaked<Eval3>(baked, fss3, natural, opts.reps));
        printRow("D relperm baked", "eval6", runRelpermBaked<Eval6>(baked, fss6, natural, opts.reps));
        printRow("D relperm baked", "shuffled3", runRelpermBaked<Eval3>(baked, fss3, shuffled, opts.reps));
        {
            const auto st = baked.storage();
            std::printf("storage: pool %zu B (%.1f KB, %zu tables, random access) | "
                        "per-cell %zu B (%.2f MB, %zu cells x %zu B, streamed) | "
                        "distinct factor sets %zu\n",
                        st.poolBytes, st.poolBytes / 1024.0, baked.numTables(),
                        st.perCellBytes, st.perCellBytes / 1048576.0,
                        st.nc, st.nc ? st.perCellBytes / st.nc : 0,
                        st.distinctFactorSets);
        }
        std::printf("uniform re-tabulation: %zu nodes (%.1f KB), max dev %.2e\n",
                    baked.numUniformNodes(), baked.numUniformNodes() * 4 * 8 / 1024.0,
                    baked.maxUniformError());
        printRow("E relperm uniform (no search)", "double",
                 runRelpermBaked<Scalar, true>(baked, fssD, natural, opts.reps));
        printRow("E relperm uniform (no search)", "eval3",
                 runRelpermBaked<Eval3, true>(baked, fss3, natural, opts.reps));
        printRow("E relperm uniform (no search)", "shuffled3",
                 runRelpermBaked<Eval3, true>(baked, fss3, shuffled, opts.reps));
    }
    else {
        std::printf("D relperm: skipped (multiplexer approach not Default)\n");
    }
    if (concretePvtOk) {
        constexpr auto WG = Opm::GasPvtApproach::WetGas;
        constexpr auto DG = Opm::GasPvtApproach::DryGas;
        printRow("D pvt analytic", "double",
                 wetGas ? runPvtAnalytic<Scalar, WG>(fssD, pvtGroups, opts.reps)
                        : runPvtAnalytic<Scalar, DG>(fssD, pvtGroups, opts.reps));
        printRow("D pvt analytic", "eval3",
                 wetGas ? runPvtAnalytic<Eval3, WG>(fss3, pvtGroups, opts.reps)
                        : runPvtAnalytic<Eval3, DG>(fss3, pvtGroups, opts.reps));
        printRow("D pvt analytic", "eval6",
                 wetGas ? runPvtAnalytic<Eval6, WG>(fss6, pvtGroups, opts.reps)
                        : runPvtAnalytic<Eval6, DG>(fss6, pvtGroups, opts.reps));
    }
    else {
        std::printf("D pvt: skipped (needs LiveOil+Wet/DryGas+CCWater)\n");
    }

    // --- relperm with hysteresis enabled (master-comparable) ---------------
    // Everything above runs hysteresis-free traits, as the experimental
    // drivers do. Master honours EHYSTR, which inserts EclHysteresisTwoPhaseLaw
    // into the stack; this section measures the same two variants with that
    // law present, at its initial (drainage) state.
    {
        using HystTraits = Opm::ThreePhaseMaterialTraits<Scalar,
                                                         FluidSystem::waterPhaseIdx,
                                                         FluidSystem::oilPhaseIdx,
                                                         FluidSystem::gasPhaseIdx,
                                                         /*enableHysteresis=*/true,
                                                         /*enableEndpointScaling=*/true>;
        using HystManager = Opm::EclMaterialLaw::Manager<HystTraits>;
        using HystLaw = typename HystManager::MaterialLaw;

        HystManager hystManager;
        hystManager.initFromState(eclState);
        {
            std::function<std::vector<int>(const Opm::FieldPropsManager&, const std::string&, bool)>
                lookup2 = [](const Opm::FieldPropsManager& fp, const std::string& prop, bool translate) {
                    std::vector<int> dest;
                    const auto& raw = fp.get_int(prop);
                    dest.resize(raw.size());
                    for (std::size_t i = 0; i < raw.size(); ++i) {
                        dest[i] = raw[i] - static_cast<int>(translate);
                    }
                    return dest;
                };
            std::function<unsigned(unsigned)> id2 = [](unsigned i) { return i; };
            hystManager.initParamsForElements(eclState, eclState.getInputGrid().getNumActive(),
                                              lookup2, id2);
        }
        std::printf("\n== relperm with hysteresis traits (enableHysteresis=%d) ==\n",
                    int(hystManager.enableHysteresis()));
        printRow("generic + hysteresis law", "eval3",
                 runRelperm<Eval3, HystLaw, HystManager>(hystManager, fss3, natural, opts.reps));
        // Same arithmetic, one shared (L1-resident) params object: isolates the
        // per-cell parameter footprint from the dispatch/arithmetic cost.
        printRow("generic + hysteresis, pinned params", "eval3",
                 runRelperm<Eval3, HystLaw, HystManager>(hystManager, fss3, natural, opts.reps, true));
        printRow("flattened + hysteresis law", "eval3",
                 runRelpermFlattened<Eval3, HystLaw, HystManager>(hystManager, fss3, natural, opts.reps));

        using HystBaked = Opm::Bench::BakedRelperm<HystLaw, HystManager, Scalar>;
        HystBaked hb;
        hb.build(hystManager, nc);
        const auto hs = hb.storage();
        std::printf("baked+hyst: %zu tables, pool %.1f KB, per-cell %zu B, max bake err %.2e\n",
                    hb.numTables(), hs.poolBytes / 1024.0,
                    hs.nc ? hs.perCellBytes / hs.nc : 0, hb.maxBakeError());
        printRow("baked + hysteresis (drainage)", "eval3",
                 runRelpermBakedT<Eval3, HystBaked>(hb, fss3, natural, opts.reps));
        // Baseline for the imbibition check below: any disagreement already
        // present here is not a hysteresis problem.
        std::printf("drainage-state ");
        compareBakedVsGeneric<Eval3, HystLaw, HystManager, HystBaked>(hystManager, hb, fss3, nc);
        // force every cell onto the imbibition branch to price the reversed path
        hb.forceReversedForBenchmark();
        printRow("baked + hysteresis (all reversed)", "eval3",
                 runRelpermBakedT<Eval3, HystBaked>(hb, fss3, natural, opts.reps));

        // Genuine reversal: replay a saturation history through the manager's
        // own updateHysteresis so both two-phase systems record a turning
        // point, then evaluate at the original states. Cells now take the
        // imbibition branch for real, and baked must still match the law.
        driveIntoImbibition<HystManager>(hystManager, cells, nc);
        hb.refreshHysteresisState(hystManager, nc);
        compareBakedVsGeneric<Eval3, HystLaw, HystManager, HystBaked>(hystManager, hb, fss3, nc);
        printRow("generic + hysteresis (imbibition)", "eval3",
                 runRelperm<Eval3, HystLaw, HystManager>(hystManager, fss3, natural, opts.reps));
        printRow("baked + hysteresis (imbibition)", "eval3",
                 runRelpermBakedT<Eval3, HystBaked>(hb, fss3, natural, opts.reps));
    }

    // --- PVT composition+b/mu: recomputed vs shared segment index ----------
    if (concretePvtOk) {
        constexpr auto WG = Opm::GasPvtApproach::WetGas;
        constexpr auto DG = Opm::GasPvtApproach::DryGas;
        std::printf("\n== pvt composition+b/mu (IQ sequence), eval3 ==\n");
        printRow("recomputed RsSat (master)", "eval3",
                 wetGas ? runPvtIQLike<Eval3, WG, false>(fss3, pvtGroups, opts.reps)
                        : runPvtIQLike<Eval3, DG, false>(fss3, pvtGroups, opts.reps));
        printRow("shared segment index", "eval3",
                 wetGas ? runPvtIQLike<Eval3, WG, true>(fss3, pvtGroups, opts.reps)
                        : runPvtIQLike<Eval3, DG, true>(fss3, pvtGroups, opts.reps));
        printRow("shared + skip when no gas", "eval3",
                 wetGas ? runPvtIQLike<Eval3, WG, true, true>(fss3, pvtGroups, opts.reps)
                        : runPvtIQLike<Eval3, DG, true, true>(fss3, pvtGroups, opts.reps));
    }

    // --- mixed-loop full IQ with baked relperm ------------------------------
    // Same fast IQ, same single per-cell loop as the simulator, but the
    // relperm block goes through the baked tables (compile-time hook). This is
    // the missing cell of the (mixed/separate x generic/baked) matrix.
    if (defaultMultiplexer && nornePvtCombo) {
        Baked bakedFB;
        bakedFB.build(manager, nc, /*force=*/true);
        using FastBakedTag = Opm::Properties::TTag::BenchIQFastBaked;
        StubIQProblemBaked<FastBakedTag> problemFB(manager, poro, rockComp, rockRefP,
                                                   benchTemperature);
        problemFB.setBaked(&bakedFB);
        const auto pvsFB = makePrimaryVariables<FastBakedTag>(states, pvtnum);
        std::printf("\n== full IQ, mixed loop, baked relperm inside ==\n");
        printRow("B+baked fast IQ", "natural",
                 runIQVariant<FastBakedTag>(problemFB, pvsFB, natural, opts.reps));
        printRow("B+baked fast IQ", "sorted",
                 runIQVariant<FastBakedTag>(problemFB, pvsFB, sorted, opts.reps));
        printRow("B+baked fast IQ", "shuffled",
                 runIQVariant<FastBakedTag>(problemFB, pvsFB, shuffled, opts.reps));
    }

    return 0;
}
