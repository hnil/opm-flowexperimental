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
#include <opm/models/utils/start.hh>
#include <opm/simulators/flow/FlowProblem.hpp>
#include "eclnewtonmethod.hh"
#include "ebos.hh"
#include <opm/simulators/flow/Main.hpp>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/flowexperimental/blackoilintensivequantitiessimple.hh>
#include "BlackOilModelFvNoCache.hpp"

// #include <opm/models/common/diffusionmodule.hh>
#include <tests/problems/co2ptflashproblem.hh>

// // the current code use eclnewtonmethod adding other conditions to proceed_ should do the trick for KA
// // adding linearshe sould be chaning the update_ function in the same class with condition that the error is reduced.
// the trick is to be able to recalculate the residual from here.
// unsure where the timestepping is done from suggestedtime??
// suggestTimeStep is taken from newton solver in problem.limitTimestep
namespace Opm{
    template<typename TypeTag>
    class OutputAuxModule : public BaseAuxiliaryModule<TypeTag>
    {

    };

}


namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowProblemEbos {
    using InheritsFrom = std::tuple<EbosTypeTag>;
};
}

template<class TypeTag>
struct Model<TypeTag, TTag::EclFlowProblemEbos> {
    using type = BlackOilModelFvNoCache<TypeTag>;
};
template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemEbos> {
     using type = BlackOilIntensiveQuantitiesSimple<TypeTag>;
};
// Set the problem class
template<class TypeTag>
struct Problem<TypeTag, TTag::EclFlowProblemEbos> {
    using type = EbosProblem<TypeTag>;
};


template<class TypeTag>
struct ThreadsPerProcess<TypeTag, TTag::EclFlowProblemEbos> {
    static constexpr int value = 1;
};

template<class TypeTag>
struct ContinueOnConvergenceError<TypeTag, TTag::EclFlowProblemEbos> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct EclNewtonSumTolerance<TypeTag, TTag::EclFlowProblemEbos> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-5;
};

// the default for the allowed volumetric error for oil per second
template<class TypeTag>
struct NewtonTolerance<TypeTag, TTag::EclFlowProblemEbos> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-2;
};

// set fraction of the pore volume where the volumetric residual may be violated during
// strict Newton iterations
template<class TypeTag>
struct EclNewtonRelaxedVolumeFraction<TypeTag, TTag::EclFlowProblemEbos> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.0;
};

template<class TypeTag>
struct EclNewtonRelaxedTolerance<TypeTag, TTag::EclFlowProblemEbos> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 10*getPropValue<TypeTag, Properties::NewtonTolerance>();
};

//template<class TypeTag>
//struct Linearizer<TypeTag, TTag::EclFlowProblemEbos> { using type = TpfaLinearizer<TypeTag>; };

// template<class TypeTag>
// struct LocalResidual<TypeTag, TTag::EclFlowProblemEbos> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::EclFlowProblemEbos> { static constexpr bool value = false; };


template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::EclFlowProblemEbos> { static constexpr bool value = false; };

//static constexpr bool has_disgas_in_water = getPropValue<TypeTag, Properties::EnableDisgasInWater>();

template<class TypeTag>
struct Simulator<TypeTag, TTag::EclFlowProblemEbos> { using type = Opm::Simulator<TypeTag>; };

// template<class TypeTag>
// struct LinearSolverBackend<TypeTag, TTag::EclFlowProblemEbos> {
//     using type = ISTLSolver<TypeTag>;
// };

// // Set the problem class
// template<class TypeTag>
// struct Problem<TypeTag, TTag::EbosTypeTag> {
//     using type = EbosProblem<TypeTag>;
// };

// template<class TypeTag>
// struct EclEnableAquifers<TypeTag, TTag::EclFlowProblemTest> {
//     static constexpr bool value = false;
// };
}
}

// CO2
namespace Opm::Properties {
    // namespace TTag {
    //      struct CO2PTBaseProblem {};
    // }
    // end namespace TTag

   namespace TTag {
   struct CO2PTEcfvProblem {
        //using InheritsFrom = std::tuple<CO2PTBaseProblem, FlashModel>;
       //missing OutPutBlackOil
       using InheritsFrom = std::tuple<FlashModel, FlowModelParameters, VtkTracer, CpGridVanguard, FlashModel, EclTimeSteppingParameters>;
       //using InheritsFrom = std::tuple<VtkTracer, OutputBlackOil, CpGridVanguard>;
   };
   }

// template <class TypeTag>
// struct NumComp<TypeTag, TTag::CO2PTEcfvProblem> {
//     static constexpr int value = 3;
// };
#if 1
template <class TypeTag>
struct MaterialLaw<TypeTag, TTag::CO2PTEcfvProblem> {
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Traits = Opm::TwoPhaseMaterialTraits<Scalar,
                                               //  /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx, // TODO
                                               /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                               /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx>;

    // define the material law which is parameterized by effective saturation
    using EffMaterialLaw = Opm::NullMaterial<Traits>;
    //using EffMaterialLaw = Opm::BrooksCorey<Traits>;

public:
    using type = EffMaterialLaw;
    //using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;
    // //using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;

    //using type = typename EclMaterialLawManager::MaterialLaw;
};
#endif
#if 0
    template<class TypeTag>
    struct MaterialLaw<TypeTag, TTag::CO2PTEcfvProblem>
    {
    private:
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

        using Traits = ThreePhaseMaterialTraits<Scalar,
                                                /*wettingPhaseIdx=*/ -1 , //FluidSystem::waterPhaseIdx,
                                                /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                                /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx>;

    public:
        using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;
        //using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;

        using type = typename EclMaterialLawManager::MaterialLaw;
    };
    template<class TypeTag>
    struct Indices<TypeTag, TTag::CO2PTEcfvProblem>
    {
    private:
        // it is unfortunately not possible to simply use 'TypeTag' here because this leads
        // to cyclic definitions of some properties. if this happens the compiler error
        // messages unfortunately are *really* confusing and not really helpful.
        using BaseTypeTag = TTag::FlowProblem;
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


#endif

template <class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::CO2PTEcfvProblem>
{
    using type = TTag::EcfvDiscretization;
};
// TODO: the following should be from the DATA file later.
template <class TypeTag>
struct EpisodeLength<TypeTag, TTag::CO2PTEcfvProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.1 * 60. * 60.;
};

template<class TypeTag>
struct CellsX<TypeTag, TTag::CO2PTEcfvProblem> { static constexpr int value = 30; };
template<class TypeTag>
struct CellsY<TypeTag, TTag::CO2PTEcfvProblem> { static constexpr int value = 1; };
// CellsZ is not needed, while to keep structuredgridvanguard.hh compile
template<class TypeTag>
struct CellsZ<TypeTag, TTag::CO2PTEcfvProblem> { static constexpr int value = 1; };

template <class TypeTag>
struct Initialpressure<TypeTag, TTag::CO2PTEcfvProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 75.e5;
};
// TODO: the above values should be from the DATA file.

template <class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::CO2PTEcfvProblem>
{
    using type = TTag::AutoDiffLocalLinearizer;
};

// Set the problem property
template <class TypeTag>
struct Problem<TypeTag, TTag::CO2PTEcfvProblem>
{
    using type = Opm::CO2PTProblem<TypeTag>;
    //using type = Opm::FlowProblem<TypeTag>;
};

template <class TypeTag>
struct EnableBrine<TypeTag, TTag::CO2PTEcfvProblem> {
    static constexpr bool value = false;
};

template <class TypeTag>
struct FlashSolver<TypeTag, TTag::CO2PTEcfvProblem> {
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

public:
    using type = Opm::PTFlash<Scalar, FluidSystem>;
};


template <class TypeTag>
struct NumComp<TypeTag, TTag::CO2PTEcfvProblem> {
    static constexpr int value = 3;
};

// set the defaults for the problem specific properties
 template <class TypeTag>
 struct Temperature<TypeTag, TTag::CO2PTEcfvProblem> {
     using type = GetPropType<TypeTag, Scalar>;
     static constexpr type value = 423.25;//TODO
 };

template <class TypeTag>
struct SimulationName<TypeTag, TTag::CO2PTEcfvProblem> {
    static constexpr auto value = "co2_ptflash";
};


template <class TypeTag>
struct FluidSystem<TypeTag, TTag::CO2PTEcfvProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr int num_comp = getPropValue<TypeTag, Properties::NumComp>();

public:
    using type = Opm::GenericOilGasFluidSystem<Scalar, num_comp>;
};


template<class TypeTag>
struct EnableMech<TypeTag, TTag::CO2PTEcfvProblem> {
    static constexpr bool value = false;
};
// disable all extensions supported by black oil model. this should not really be
// necessary but it makes things a bit more explicit
template<class TypeTag>
struct EnablePolymer<TypeTag, TTag::CO2PTEcfvProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableSolvent<TypeTag, TTag::CO2PTEcfvProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::CO2PTEcfvProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableFoam<TypeTag, TTag::CO2PTEcfvProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableExtbo<TypeTag, TTag::CO2PTEcfvProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableMICP<TypeTag, TTag::CO2PTEcfvProblem> {
    static constexpr bool value = false;
};

// disable thermal flux boundaries by default
template<class TypeTag>
struct EnableThermalFluxBoundaries<TypeTag, TTag::CO2PTEcfvProblem> {
    static constexpr bool value = false;
};


template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::CO2PTEcfvProblem>
{
    static constexpr bool value = false;
};


} // namespace Opm::Properties





int main(int argc, char** argv)
{
    //using TypeTag = Opm::Properties::TTag::EclFlowProblemEbos;
    using TypeTag = Opm::Properties::TTag::CO2PTEcfvProblem;
    Opm::registerEclTimeSteppingParameters<TypeTag>();
    return Opm::start<TypeTag>(argc, argv);
}
