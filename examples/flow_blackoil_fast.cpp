/*
  Copyright 2020, NORCE AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "config.h"
#include "tracy/Tracy.hpp"
#include "tracy/TracyC.h"
#define OPM_TIMEBLOCK(blockname) ZoneNamedN(blockname, #blockname, true);
#define OPM_TIMEBLOCK_LOCAL(blockname);// ZoneNamedN(blockname, #blockname, true);
#include <opm/simulators/flow/Main.hpp>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/flowexperimental/eclblackoilintensivequantities.hh>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManagerTable.hpp>

// initialization modifications to be able to inititialize with new material manager
#include <ebos/equil/equilibrationhelpers.hh>
#include <ebos/equil/equilibrationhelpers_impl.hh>//new file in flowexperimental
#include <ebos/equil/initstateequil.hh>
#include <ebos/equil/initstateequil_impl.hh>//new file in flow experimental

namespace Opm{
    template<typename TypeTag>
    class EclProblemNew: public EclProblem<TypeTag>{
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using DirectionalMobilityPtr = ::Opm::Utility::CopyablePtr<DirectionalMobility<TypeTag, Evaluation>>;
        using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        enum { numPhases = FluidSystem::numPhases };
    public:
        EclProblemNew(Simulator& simulator): EclProblem<TypeTag>(simulator){
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


namespace Opm {
namespace Properties {
namespace TTag {
    struct EclFlowProblemTest {
        using InheritsFrom = std::tuple<EclFlowProblem>;
    };
}
    template<class TypeTag>
    struct Linearizer<TypeTag, TTag::EclFlowProblemTest> { using type = TpfaLinearizer<TypeTag>; };
    
    template<class TypeTag>
    struct LocalResidual<TypeTag, TTag::EclFlowProblemTest> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

    template<class TypeTag>
    struct EnableDiffusion<TypeTag, TTag::EclFlowProblemTest> { static constexpr bool value = false; };
    template<class TypeTag>
    struct Model<TypeTag, TTag::EclFlowProblemTest> {
        using type = BlackOilModelFv<TypeTag>;
    };
    template<class TypeTag>
    struct Problem<TypeTag, TTag::EclFlowProblemTest> {
        using type = EclProblemNew<TypeTag>;
    };

    template<class TypeTag>
    struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemTest> {
    using type = EclBlackOilIntensiveQuantities<TypeTag>;
    };
    
    
    template<class TypeTag>
    struct MaterialLaw<TypeTag, TTag::EclFlowProblemTest>
    {
    private:
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        
        using Traits = ThreePhaseMaterialTraits<Scalar,
                                                /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                                /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                                /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx>;
        
    public:
        using EclMaterialLawManager = ::Opm::EclMaterialLawManagerTable<Traits>;
        //using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;
        
        using type = typename EclMaterialLawManager::MaterialLaw;
    };

};    
   
}

int main(int argc, char** argv)
{
    OPM_TIMEBLOCK(fullSimulation);
    using TypeTag = Opm::Properties::TTag::EclFlowProblemTest;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
//    return Opm::start<TypeTag>(argc, argv);
}
