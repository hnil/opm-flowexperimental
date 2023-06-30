#ifndef BLACK_OIL_MODEL_FV_HPP
#define BLACK_OIL_MODEL_FV_HPP
#include <ebos/FIBlackOilModel.hpp>
namespace Opm{
    template<typename TypeTag>
    class BlackOilModelFv: public BlackOilModel<TypeTag>{
        using Parent = BlackOilModel<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
        using Indices = GetPropType<TypeTag, Properties::Indices>;
        static constexpr bool waterEnabled = Indices::waterEnabled;
        static constexpr bool gasEnabled = Indices::gasEnabled;
        static constexpr bool oilEnabled = Indices::oilEnabled;
        enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };



    public:
        BlackOilModelFv(Simulator& simulator): BlackOilModel<TypeTag>(simulator){
        }
        void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx){
            OPM_TIMEBLOCK_LOCAL(updateIntensiveQuantities);
//            std::cout << "----------------------Update quantities-------------------\n"
//                      << std::flush;
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
                //intquant.update(problem,priVars, globalSpaceIdx, timeIdx, waterpvt, gaspvt, oilpvt);
            }
            
            std::fill(this->intensiveQuantityCacheUpToDate_[timeIdx].begin(),
                      this->intensiveQuantityCacheUpToDate_[timeIdx].end(),
                      /*value=*/true);
            
        }

        const IntensiveQuantities& intensiveQuantities(unsigned globalIdx, unsigned timeIdx) const{
            OPM_TIMEBLOCK_LOCAL(intensiveQuantities);
            const auto& primaryVars = this->solution(timeIdx);
            const auto& problem = this->simulator_.problem();
            const auto intquant = this->cachedIntensiveQuantities(globalIdx, timeIdx);
            if (!this->enableIntensiveQuantityCache_){
                OPM_THROW(std::logic_error, "Run without intentive quantites not enabled: Use --enable-intensive-quantity=true");
            }
            if(!intquant){
                OPM_THROW(std::logic_error, "Intensive quantites need to be updated in code");
            }    
            return *intquant;    
        }

    };
}
#endif
