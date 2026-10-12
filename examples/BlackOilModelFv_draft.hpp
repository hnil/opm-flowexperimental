#ifndef BLACK_OIL_MODEL_FV_HPP
#define BLACK_OIL_MODEL_FV_HPP
#include <ebos/FIBlackOilModel.hh>
namespace Opm{
    template<typename TypeTag>
    class BlackOilModelFv: public FIBlackOilModel<TypeTag>{
        using Parent = BlackOilModel<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    public:
        BlackOilModelFv(Simulator& simulator): BlackOilModel<TypeTag>(simulator){
        }
        void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx){
            std::cout << "----------------------Update quantities-------------------\n"
                      << std::flush;
//            Parent::invalidateAndUpdateIntensiveQuantities(timeIdx);
//             Parent::invalidateAndUpdateIntensiveQuantitiesSimple(*this,solution,/*timeIdx*/0);
            const auto& primaryVars = this->solution(timeIdx);
            const auto& problem = this->simulator_.problem();
            size_t numGridDof = primaryVars.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif                
            for (unsigned dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
                const auto& primaryVar = primaryVars[dofIdx];
                auto& intquant = this->intensiveQuantityCache_[timeIdx][dofIdx];
                intquant.update(problem, primaryVar, dofIdx, timeIdx);
            }
            
            std::fill(this->intensiveQuantityCacheUpToDate_[timeIdx].begin(),
                      this->intensiveQuantityCacheUpToDate_[timeIdx].end(),
                      /*value=*/true);
            
        }

        // standard flow
        const IntensiveQuantities& intensiveQuantities(unsigned globalIdx, unsigned timeIdx){
            const auto& primaryVars = this->solution(timeIdx);
            const auto& problem = this->simulator_.problem();
            auto* intquant = this->cachedIntensiveQuantities(unsigned globalIdx, unsigned timeIdx);
            if (!enableIntensiveQuantityCache_){
                OPM_THROW(std::logic_error, "Run without intentive quantites not enabled: Use --enable-intensive-quantity=true");
            }
            if(!intquant){
                OPM_THROW(std::logic_error, "Intensive quantites need to be updated in code");
            }    
            return &intquant;    
        }

        // flow with enable IntenisveQuantities future
        const IntensiveQuantities& intensiveQuantities(unsigned globalIdx, unsigned timeIdx){
            const auto& primaryVars = this->solution(timeIdx);
            const auto& problem = this->simulator_.problem();
            if (!enableIntensiveQuantityCache_){
                OPM_THROW(std::logic_error, "Run without intentive quantites not enabled: Use --enable-intensive-quantity=true");
            }else{
                const auto& primaryVars = this->solution(timeIdx);
                const auto& problem = this->simulator_.problem();
                auto& intquant = this->intensiveQuantityCache_[timeIdx][dofIdx];
                if(!intensiveQuantityCacheUpToDate_[timeIdx][globalIdx]){
                    const auto& primaryVar = primaryVars[dofIdx];
                    intquant.update(problem, primaryVar, dofIdx, timeIdx);
                    intensiveQuantityCacheUpToDate_[timeIdx][globalIdx] = true;
                }
                return intquant;
            }
        }

        // flow with intensive quantities off future
        IntensiveQuantities intensiveQuantities(unsigned globalIdx, unsigned timeIdx){
            if (enableIntensiveQuantityCache_){
                const auto& primaryVars = this->solution(timeIdx);
                const auto& problem = this->simulator_.problem();
                auto& intquant = this->intensiveQuantityCache_[timeIdx][dofIdx];
                if(!intensiveQuantityCacheUpToDate_[timeIdx][globalIdx]){
                    const auto& primaryVar = primaryVars[dofIdx];
                    intquant.update(problem, primaryVar, dofIdx, timeIdx);
                }
                return intquant;
            }else{
                const auto& primaryVars = this->solution(timeIdx);
                const auto& problem = this->simulator_.problem();
                IntensiveQuantites inquant;
                const auto& primaryVar = primaryVars[dofIdx];
                intquant.update(problem, primaryVar, dofIdx, timeIdx);
            return intquant;            
        }
       
    };
}
#endif
