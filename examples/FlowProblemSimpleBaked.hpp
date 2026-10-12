// FlowProblemSimpleFast + baked relperm tables: at finishInit() the per-cell
// material law (EPS included) is baked into deduplicated piecewise-linear
// tables (see material/BakedRelperm.hpp); the fast IQ then evaluates relperm
// and capillary pressure through evalBakedRelperm() instead of the law stack.
#ifndef FLOW_PROBLEM_SIMPLE_BAKED_HPP
#define FLOW_PROBLEM_SIMPLE_BAKED_HPP

#include "FlowProblemSimpleFast.hpp"
#include "material/BakedRelperm.hpp"

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

namespace Opm {

template <typename TypeTag>
class FlowProblemSimpleBaked : public FlowProblemSimpleFast<TypeTag>
{
    using Parent = FlowProblemSimpleFast<TypeTag>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Manager = typename GetProp<TypeTag, Properties::MaterialLaw>::EclMaterialLawManager;

    enum { numPhases = FluidSystem::numPhases };
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr unsigned oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;

public:
    explicit FlowProblemSimpleBaked(Simulator& simulator)
        : Parent(simulator)
    {}

    void finishInit()
    {
        Parent::finishInit();
        const auto& manager = *this->materialLawManager();
        if (manager.hasDirectionalRelperms() || manager.hasDirectionalImbnum()) {
            throw std::runtime_error("baked relperm: directional relperms not supported");
        }
        // Deck-requested hysteresis is already ignored by this driver's
        // hysteresis-free Traits; the bake matches that law exactly. Only a
        // law that actually evolves its params (Traits hysteresis) is unsupported.
        static_assert(!MaterialLaw::Traits::enableHysteresis,
                      "baked relperm: hysteresis-enabled traits not supported yet");
        const std::size_t nc = this->simulator().model().numGridDof();
        baked_.build(manager, nc);
        std::printf("baked relperm: %zu cells -> %zu unique tables, %zu nodes (%.1f KB), "
                    "max bake err %.2e\n",
                    nc, baked_.numTables(), baked_.numNodes(),
                    baked_.numNodes() * 5 * 8 / 1024.0, baked_.maxBakeError());
        if (baked_.maxBakeError() > 1.0e-8) {
            throw std::runtime_error("baked relperm: bake validation failed, max rel err "
                                     + std::to_string(baked_.maxBakeError()));
        }
        // same-binary A/B switch: OPM_USE_BAKED=0 falls back to the generic law
        const char* env = std::getenv("OPM_USE_BAKED");
        useBaked_ = !(env != nullptr && env[0] == '0');
        // per-call hook timing (~20 ns/call overhead): OPM_HOOK_TIMING=1
        const char* timing = std::getenv("OPM_HOOK_TIMING");
        hookTiming_ = (timing != nullptr && timing[0] == '1');
        const char* fastExp = std::getenv("OPM_FAST_EXPLICIT");
        useFastExplicit_ = !(fastExp != nullptr && fastExp[0] == '0');
    }

    ~FlowProblemSimpleBaked()
    {
        const long long calls = hookCalls_.load();
        if (calls > 0) {
            std::printf("relperm hook (%s): %lld calls, %.3f s total, %.1f ns/call\n",
                        useBaked_ ? "baked" : "generic", calls,
                        hookNs_.load() * 1.0e-9, double(hookNs_.load()) / double(calls));
        }
    }

    // A/B switch for the experimental explicit-quantities rewrite that
    // FlowProblemSimpleFast installs (it drives the pre/post-step bucket):
    // OPM_FAST_EXPLICIT=0 runs master's implementation instead.
    void updateExplicitQuantities_(int episodeIdx, int timeStepSize,
                                   const bool first_step_after_restart) override
    {
        if (useFastExplicit_) {
            Parent::updateExplicitQuantities_(episodeIdx, timeStepSize, first_step_after_restart);
        }
        else {
            FlowProblemBlackoil<TypeTag>::updateExplicitQuantities_(first_step_after_restart);
        }
    }

    template <class Ev>
    void evalBakedRelperm(unsigned globalSpaceIdx,
                          const Ev& sw,
                          const Ev& sg,
                          std::array<Ev, numPhases>& mobility,
                          std::array<Ev, numPhases>& pC) const
    {
        if (!hookTiming_) {
            evalBakedRelpermImpl(globalSpaceIdx, sw, sg, mobility, pC);
            return;
        }
        const auto t0 = std::chrono::steady_clock::now();
        evalBakedRelpermImpl(globalSpaceIdx, sw, sg, mobility, pC);
        const auto dt = std::chrono::steady_clock::now() - t0;
        hookNs_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(dt).count(),
                          std::memory_order_relaxed);
        hookCalls_.fetch_add(1, std::memory_order_relaxed);
    }

    template <class Ev>
    void evalBakedRelpermImpl(unsigned globalSpaceIdx,
                              const Ev& sw,
                              const Ev& sg,
                              std::array<Ev, numPhases>& mobility,
                              std::array<Ev, numPhases>& pC) const
    {
        if (!useBaked_) {
            Bench::SatOnlyFluidState<Ev, numPhases> fs;
            fs.sat[waterPhaseIdx] = sw;
            fs.sat[gasPhaseIdx] = sg;
            fs.sat[oilPhaseIdx] = 1.0 - sw - sg;
            const auto& materialParams = this->materialLawParams(globalSpaceIdx);
            MaterialLaw::capillaryPressures(pC, materialParams, fs);
            MaterialLaw::relativePermeabilities(mobility, materialParams, fs);
            return;
        }
        // kr and pc differences written in place; pC relative to oil is all
        // the pressure update consumes.
        baked_.evaluate(globalSpaceIdx, sw, sg, mobility,
                        pC[waterPhaseIdx], pC[gasPhaseIdx]);
        pC[oilPhaseIdx] = 0.0;
    }

    const Bench::BakedRelperm<MaterialLaw, Manager, Scalar>& bakedRelperm() const
    { return baked_; }

private:
    Bench::BakedRelperm<MaterialLaw, Manager, Scalar> baked_;
    bool useBaked_ = true;
    bool hookTiming_ = false;
    bool useFastExplicit_ = true;
    mutable std::atomic<long long> hookNs_{0};
    mutable std::atomic<long long> hookCalls_{0};
};

} // namespace Opm

#endif
