// Minimal problem stub so BlackOilIntensiveQuantities::update(problem, priVars,
// globalIdx, timeIdx) can run outside a simulator (see bench_iq_evaluation.cpp).
// Modeled on opm-simulators/tests/test_blackoilprimaryvariables.cpp, extended
// with the members the context-free IQ update actually calls.
#ifndef OPM_BENCH_STUB_IQ_PROBLEM_HPP
#define OPM_BENCH_STUB_IQ_PROBLEM_HPP

#include <opm/models/discretization/common/linearizationtype.hh>
#include <opm/models/common/directionalmobility.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/utility/CopyablePtr.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

#include <array>
#include <cassert>
#include <vector>

namespace Opm::Bench {

struct StubLinearizer
{
    LinearizationType linearizationType_{};
    const LinearizationType& getLinearizationType() const
    { return linearizationType_; }
};

struct StubModel
{
    StubLinearizer linearizer_{};
    const StubLinearizer& linearizer() const
    { return linearizer_; }
};

template <class TypeTag>
class StubIQProblem
{
public:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using DirectionalMobilityPtr = Utility::CopyablePtr<DirectionalMobility<TypeTag>>;

    static constexpr int numPhases = FluidSystem::numPhases;

    // The manager whose per-cell params back materialLawParams(). Traits must
    // match the TypeTag's MaterialLaw property.
    using Traits = typename MaterialLaw::Traits;
    using Manager = EclMaterialLaw::Manager<Traits>;

    StubIQProblem(const Manager& manager,
                  std::vector<Scalar> porosity,
                  Scalar rockCompressibility,
                  Scalar rockReferencePressure,
                  Scalar temperature)
        : manager_(manager)
        , porosity_(std::move(porosity))
        , rockComp_(rockCompressibility)
        , rockRefP_(rockReferencePressure)
        , temperature_(temperature)
    {}

    const StubModel& model() const
    { return model_; }

    unsigned satnumRegionIndex(unsigned elemIdx) const
    { return manager_.satnumRegionIdx(elemIdx); }

    Scalar temperature(unsigned, unsigned) const
    { return temperature_; }

    const MaterialLawParams& materialLawParams(unsigned elemIdx) const
    { return manager_.materialLawParams(elemIdx); }

    Scalar porosity(unsigned elemIdx, unsigned) const
    { return porosity_[elemIdx]; }

    Scalar rockCompressibility(unsigned) const
    { return rockComp_; }

    Scalar rockReferencePressure(unsigned) const
    { return rockRefP_; }

    template <class Ev, class IQ>
    Ev rockCompPoroMultiplier(const IQ&, unsigned) const
    { return 1.0; }

    template <class Ev, class IQ>
    Ev rockCompTransMultiplier(const IQ&, unsigned) const
    { return 1.0; }

    // Sentinels: never limit Rs/Rv, so the saturated tables decide.
    Scalar maxGasDissolutionFactor(unsigned, unsigned) const
    { return 1.0e30; }

    Scalar maxOilVaporizationFactor(unsigned, unsigned) const
    { return 1.0e30; }

    Scalar maxOilSaturation(unsigned) const
    { return 1.0; }

    // Mirrors FlowProblem::updateRelperms (FlowProblem.hpp:860), minus the
    // directional branch (asserted off at harness init).
    template <class FluidState, class ...Args>
    void updateRelperms(std::array<Evaluation, numPhases>& mobility,
                        DirectionalMobilityPtr&,
                        FluidState& fluidState,
                        unsigned globalSpaceIdx) const
    {
        using ContainerT = std::array<Evaluation, numPhases>;
        const auto& materialParams = materialLawParams(globalSpaceIdx);
        MaterialLaw::template relativePermeabilities<ContainerT, FluidState, Args...>
            (mobility, materialParams, fluidState);
    }

private:
    const Manager& manager_;
    std::vector<Scalar> porosity_;
    Scalar rockComp_;
    Scalar rockRefP_;
    Scalar temperature_;
    StubModel model_{};
};

} // namespace Opm::Bench

#endif
