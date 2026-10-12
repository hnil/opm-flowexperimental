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
        
    protected:
        void updateExpliciteQuantities_(){
            OPM_TIMEBLOCK(updateExpliciteQuantitiesSimple);            
            const auto& primaryVars = this->modle().solution(/*timeIdx*/0);
            const auto& problem = this->simulator().problem();
            const auto& model = this->simulato().model();
            const auto& vanguard = this->simulator().vangaurd();
            size_t numGridDof = primaryVars.size();
            std::string failureMsg("Error updating explicit quantities");
            OPM_BEGIN_PARALLEL_TRY_CATCH();                       
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (unsigned dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
                bool changed=false;
                const auto& iq = *model.cachedIntensiveQuantities(dofIdx, /*timeIdx=*/ 0);
                if (this->maxWaterSaturation_.empty()){
                    const bool invalidateFromMaxWaterSat = updateMaxWaterSaturation_(dofIdx,iq);
                    changed = changed || invalidateFromMaxWaterSat;
                }
                if (this->minOilPressure_.empty()){
                    const bool invalidateFromMinPressure = updateMinPressure_(dofIdx,iq);
                    changed = changed || invalidateFromMinPressure;
                }
                
                // update hysteresis and max oil saturation used in vappars
                if (!this->materialLawManager_->enableHysteresis()){
                    //need bools in opm-material
                    const bool invalidateFromHyst = this->materialLawManager()->updateHysteresis_(dofIdx,iq);
                    changed = changed || invalidateFromHyst;
                }
                int episodeIdx = this->episodeIndex();
                if (this->vapparsActive(episodeIdx)) {
                    const bool invalidateFromMaxOilSat = updateMaxOilSaturation_(dofIdx,iq);
                    changed = changed || invalidateFromMaxOilSat;
                }
                
                // the derivatives may have change
                
                if(changed){
                    OPM_TIMEBLOCK_LOCAL(updateIntensiveQuantityHyst);
                    const auto& primaryVar = primaryVars[dofIdx];
                    //auto& intquant = this->intensiveQuantityCache_[timeIdx][dofIdx];
                    iq.update(problem, primaryVar, dofIdx,/*timeIdx*/0);
                }
                if constexpr (getPropValue<TypeTag, Properties::EnablePolymer>()){
                    updateMaxPolymerAdsorption_(dofIdx,iq);
                }
            }
            OPM_END_PARALLEL_TRY_CATCH(failureMsg, vanguard.grid().comm());
        }


        // }
//         template<class UpdateFunc>
//         void updateProperty_(const std::string& failureMsg,
//                              UpdateFunc func)
//         {
//             OPM_TIMEBLOCK(updateProperty);            
//             const auto& primaryVars = this->solution(/*timeIdx*/0);
//             const auto& problem = this->simulator_.problem();
//             const auto& model = this->simulator_.model();
//             const auto& vanguard = this->simulator_.vangaurd();
//             size_t numGridDof = primaryVars.size();
// //            OPM_BEGIN_PARALLEL_TRY_CATCH();                       
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
//             for (unsigned dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
//                 const auto& iq = *model.cachedIntensiveQuantities(dofIdx, /*timeIdx=*/ 0);
//                 func(dofIdx, iq);
//             }
// //            OPM_END_PARALLEL_TRY_CATCH(failureMsg, vanguard.grid().comm());
//         }
            
    };

}


#endif
