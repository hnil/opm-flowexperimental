#ifndef BLACK_OIL_MODEL_FV_HPP
#define BLACK_OIL_MODEL_FV_HPP

namespace Opm{
    template<typename TypeTag>
    class BlackOilModelFv: public BlackOilModel<TypeTag>{
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
    };
}
#endif
