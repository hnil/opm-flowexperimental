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
// #if HAVE_TRACY
// #include "tracy/Tracy.hpp"
// #include "tracy/TracyC.h"
// #define OPM_TIMEBLOCK(blockname);// ZoneNamedN(blockname, #blockname, true);
// #define OPM_TIMEBLOCK_LOCAL(blockname);// ZoneNamedN(blockname, #blockname, true);
// #endif
#include <opm/simulators/flow/Main.hpp>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/flowexperimental/blackoilintensivequantitiesdrygas.hh>
#include <opm/flowexperimental/blackoilintensivequantitiesdrygas.hh>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManagerTable.hpp>
#include "BlackOilModelFvFast.hpp"

namespace Opm {
namespace Properties {
namespace TTag {
    struct EclFlowProblemTest {
        using InheritsFrom = std::tuple<FlowProblem>;
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
        using type = BlackOilModelFvFast<TypeTag>;
    };

    template<class TypeTag>
    struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemTest> {
    using type = BlackOilIntensiveQuantitiesDryGas<TypeTag>;
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
        using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;

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
