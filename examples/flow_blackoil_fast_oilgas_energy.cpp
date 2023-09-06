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
#if HAVE_TRACY
#include "tracy/Tracy.hpp"
#include "tracy/TracyC.h"
#define OPM_TIMEBLOCK(blockname) ZoneNamedN(blockname, #blockname, true);
#define OPM_TIMEBLOCK_LOCAL(blockname);// ZoneNamedN(blockname, #blockname, true);
#endif
#include <opm/simulators/flow/Main.hpp>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/flowexperimental/eclblackoilintensivequantities.hh>
#include <opm/flowexperimental/blackoilintensivequantitiessimple.hh>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManagerTable.hpp>

// initialization modifications to be able to inititialize with new material manager
#include <ebos/equil/equilibrationhelpers.hh>
#include <ebos/equil/equilibrationhelpers_impl.hh>//new file in flowexperimental
#include <ebos/equil/initstateequil.hh>
#include <ebos/equil/initstateequil_impl.hh>//new file in flow experimental
#include "BlackOilModelFv.hpp"
#include "EclProblemSimple.hpp"

namespace Opm {
namespace Properties {
namespace TTag {
    struct EclFlowProblemTest {
        using InheritsFrom = std::tuple<EclFlowProblem>;
    };
}
    // template<class TypeTag>
    // struct Linearizer<TypeTag, TTag::EclFlowProblemTest> { using type = TpfaLinearizer<TypeTag>; };

    template<class TypeTag>
    struct LocalResidual<TypeTag, TTag::EclFlowProblemTest> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

    template<class TypeTag>
    struct EnableEnergy<TypeTag, TTag::EclFlowProblemTest> {
        static constexpr bool value = true;
    };

    template<class TypeTag>
    struct EnableDiffusion<TypeTag, TTag::EclFlowProblemTest> { static constexpr bool value = false; };
    template<class TypeTag>
    struct Model<TypeTag, TTag::EclFlowProblemTest> {
        using type = BlackOilModelFv<TypeTag>;
    };
    template<class TypeTag>
    struct Problem<TypeTag, TTag::EclFlowProblemTest> {
        using type = EclProblemSimple<TypeTag>;
    };

    template<class TypeTag>
    struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemTest> {
        using type = BlackOilIntensiveQuantitiesSimple<TypeTag>;
        //using type = EclBlackOilIntensiveQuantities<TypeTag>;

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
        using EclMaterialLawManager = ::Opm::EclMaterialLawManagerTable<Traits,2>;
        //using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;

        using type = typename EclMaterialLawManager::MaterialLaw;
    };
    template<class TypeTag>
    struct Indices<TypeTag, TTag::EclFlowProblemTest>
    {
    private:
        // it is unfortunately not possible to simply use 'TypeTag' here because this leads
        // to cyclic definitions of some properties. if this happens the compiler error
        // messages unfortunately are *really* confusing and not really helpful.
        using BaseTypeTag = TTag::EclFlowProblem;
        using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;

    public:
        typedef BlackOilTwoPhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                        getPropValue<TypeTag, Properties::EnableExtbo>(),
                                        getPropValue<TypeTag, Properties::EnablePolymer>(),
                                        getPropValue<TypeTag, Properties::EnableEnergy>(),
                                        getPropValue<TypeTag, Properties::EnableFoam>(),
                                        getPropValue<TypeTag, Properties::EnableBrine>(),
                                        /*PVOffset=*/0,
                                        /*disabledCompIdx=*/FluidSystem::waterCompIdx,
                                        getPropValue<TypeTag, Properties::EnableMICP>()> type;
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
