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
#include "co2ptflowproblem.hh"
#include <opm/simulators/flow/FlowGenericProblem.hpp>
#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>


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
    template<typename TypeTag>
    class EmptyModel : public BaseAuxiliaryModule<TypeTag>
    {
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
        using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    public:
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        EmptyModel(Simulator& /*simulator*/){};
        void init(){}
        template<class Something>
        void init(Something /*A*/){}
        void prepareTracerBatches(){};
        using NeighborSet = std::set<unsigned>;
        void linearize(SparseMatrixAdapter& /*matrix*/, GlobalEqVector& /*residual*/){};
        unsigned numDofs() const{return 0;};
        void addNeighbors(std::vector<NeighborSet>& /*neighbors*/) const{};
        //void applyInitial(){};
        void initialSolutionApplied(){};
        //void initFromRestart(const data::Aquifers& aquiferSoln);
        template <class Restarter>
        void serialize(Restarter& /*res*/){};

        template <class Restarter>
        void deserialize(Restarter& /*res*/){};

        void beginEpisode(){};
        void beginTimeStep(){};
        void beginIteration(){};
        // add the water rate due to aquifers to the source term.
        template<class RateVector, class Context>
        void addToSource(RateVector& rates, const Context& context, unsigned spaceIdx, unsigned timeIdx) const{};
        template<class RateVector>
        void addToSource(RateVector& rates, unsigned globalSpaceIdx, unsigned timeIdx) const{};
        void endIteration()const{};
        void endTimeStep(){};
        void endEpisode(){};
        void applyInitial(){};
        template<class RateType>
        void computeTotalRatesForDof(RateType& /*rate*/, unsigned /*globalIdx*/) const{};
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
       //using InheritsFrom = std::tuple<FlashModel,EbosTypeTag>;
       //missing OutPutBlackOil
       using InheritsFrom = std::tuple<FlashModel, FlowModelParameters, VtkTracer, CpGridVanguard, FlashModel, EclTimeSteppingParameters>;
       //using InheritsFrom = std::tuple<VtkTracer, OutputBlackOil, CpGridVanguard>;
   };
   }
    template<class TypeTag, class MyTypeTag>
    struct ExpliciteRockCompaction{
        using type = UndefinedProperty;
    };
// template <class TypeTag>
// struct NumComp<TypeTag, TTag::CO2PTEcfvProblem> {
//     static constexpr int value = 3;
// };
#if 0
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

    template<class TypeTag>
    struct MaterialLaw<TypeTag, TTag::CO2PTEcfvProblem>
    {
    private:
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using Indices = GetPropType<TypeTag, Properties::Indices>;

        using Traits = ThreePhaseMaterialTraits<Scalar,
                                                /*wettingPhaseIdx=*/0,
                                                /*nonWettingPhaseIdx=*/1,
                                                /*gasPhaseIdx=*/2>;

    public:
        using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;
        //using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;

        using type = typename EclMaterialLawManager::MaterialLaw;
    };

// Set the problem property
template <class TypeTag>
struct Problem<TypeTag, TTag::CO2PTEcfvProblem>
{
    //using type = Opm::CO2PTProblem<TypeTag>;
    using type = EbosProblem<TypeTag>;
};

template<class TypeTag>
struct AquiferModel<TypeTag, TTag::CO2PTEcfvProblem> {
    using type = EmptyModel<TypeTag>;
};

template<class TypeTag>
struct WellModel<TypeTag, TTag::CO2PTEcfvProblem> {
    using type = EmptyModel<TypeTag>;
};
template<class TypeTag>
struct TracerModelDef<TypeTag, TTag::CO2PTEcfvProblem> {
    using type = EmptyModel<TypeTag>;
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
// disable thermal flux boundaries by default
template<class TypeTag>
struct EnableThermalFluxBoundaries<TypeTag, TTag::CO2PTEcfvProblem> {
    static constexpr bool value = false;
};

template <class TypeTag>
struct EpisodeLength<TypeTag, TTag::CO2PTEcfvProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.1 * 60. * 60.;
};

template <class TypeTag>
struct Initialpressure<TypeTag, TTag::CO2PTEcfvProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 75.e5;
};

template <class TypeTag>
struct DomainSizeZ<TypeTag, TTag::CO2PTEcfvProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};

template<class TypeTag>
struct CellsX<TypeTag, TTag::CO2PTEcfvProblem> { static constexpr int value = 30; };
template<class TypeTag>
struct CellsY<TypeTag, TTag::CO2PTEcfvProblem> { static constexpr int value = 1; };
// CellsZ is not needed, while to keep structuredgridvanguard.hh compile
template<class TypeTag>
struct CellsZ<TypeTag, TTag::CO2PTEcfvProblem> { static constexpr int value = 1; };

} // namespace Opm::Properties





int main(int argc, char** argv)
{
    //using TypeTag = Opm::Properties::TTag::EclFlowProblemEbos;
    using TypeTag = Opm::Properties::TTag::CO2PTEcfvProblem;
    Opm::registerEclTimeSteppingParameters<TypeTag>();
    return Opm::start<TypeTag>(argc, argv);
}
