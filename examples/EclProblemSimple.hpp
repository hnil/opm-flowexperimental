#ifndef ECL_PROBLEM_OIL_SIMPLE_HPP
#define ECL_PROBLEM_OIL_SIMPLE_HPP
namespace Opm{
    template<typename TypeTag>
    class EclProblemSimple: public EclProblem<TypeTag>{
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using DirectionalMobilityPtr = ::Opm::Utility::CopyablePtr<DirectionalMobility<TypeTag, Evaluation>>;
        using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        enum { numPhases = FluidSystem::numPhases };
    public:
        EclProblemSimple(Simulator& simulator): EclProblem<TypeTag>(simulator){
        }
        template <class FluidState>
        void updateRelperms(
            std::array<Evaluation,numPhases> &mobility,
            DirectionalMobilityPtr &/*dirMob*/,
            FluidState &fluidState,
            unsigned globalSpaceIdx) const
        {
            OPM_TIMEBLOCK_LOCAL(updateRelperms);
            {
                // calculate relative permeabilities. note that we store the result into the
                // mobility_ class attribute. the division by the phase viscosity happens later.
                const auto& materialParams = this->materialLawParams(globalSpaceIdx);
                MaterialLaw::relativePermeabilities(mobility, materialParams, fluidState);
                Valgrind::CheckDefined(mobility);
            }
        };
    };

}


#endif
