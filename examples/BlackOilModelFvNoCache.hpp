#ifndef BLACK_OIL_MODEL_FV_NOCACHE_HPP
#define BLACK_OIL_MODEL_FV_NOCACHE_HPP
#include <ebos/FIBlackOilModel.hpp>
namespace Opm{
    template<typename TypeTag>
    class BlackOilModelFvNoCache: public BlackOilModel<TypeTag>{
        using Parent = BlackOilModel<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    public:
        BlackOilModelFvNoCache(Simulator& simulator): BlackOilModel<TypeTag>(simulator){
        }
        void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx){
            OPM_TIMEBLOCK(updateIntensiveQuantitiesNoCache);
//            std::cout << "----------------------Update quantities-------------------\n"
//                      << std::flush;
//            Parent::invalidateAndUpdateIntensiveQuantities(timeIdx);
//             Parent::invalidateAndUpdateIntensiveQuantitiesSimple(*this,solution,/*timeIdx*/0);
            if(!this->enableIntensiveQuantityCache_){
                return;
            }
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

        IntensiveQuantities intensiveQuantities(unsigned globalIdx, unsigned timeIdx) const{
            OPM_TIMEBLOCK_LOCAL(intensiveQuantitiesNoCache);
            const auto& primaryVar = this->solution(timeIdx)[globalIdx];
            const auto& problem = this->simulator_.problem();
            //IntensiveQuantities* intQuants = &(this->intensiveQuantityCache_[timeIdx][globalIdx]);
            if (!(this->enableIntensiveQuantityCache_) ||
                !(this->intensiveQuantityCacheUpToDate_[timeIdx][globalIdx])){
                IntensiveQuantities intQuants;
                intQuants.update(problem,primaryVar, globalIdx, timeIdx);
                return intQuants;// reqiored for updating extrution factor
            }else{
                IntensiveQuantities intQuants = (this->intensiveQuantityCache_[timeIdx][globalIdx]);
                return intQuants;
            }

        }

    };
}
#endif
