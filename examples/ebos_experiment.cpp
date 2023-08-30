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
#include <exception>
#include <ebos/eclproblem.hh>
#include <ebos/eclnewtonmethod.hh>
#include <ebos/ebos.hh>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/flowexperimental/blackoilintensivequantitiessimple.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include "BlackOilModelFvNoCache.hpp"
namespace Opm{
    template<typename TypeTag>
    class MonitoringAuxModule : public BaseAuxiliaryModule<TypeTag>
    {      
        using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
        using NeighborSet = typename BaseAuxiliaryModule<TypeTag>::NeighborSet;
        using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    public:
        void postSolve(GlobalEqVector& deltaX){
            std::cout << "Dummy PostSolve Aux" << std::endl;
        }
        void addNeighbors(std::vector<NeighborSet>& neighbors) const{};
        void applyInitial(){};
        unsigned numDofs() const{return 0;};
        void linearize(SparseMatrixAdapter& matrix, GlobalEqVector& residual){
            std::cout << "Dummy Linearize Aux" << std::endl;
        };
    };

    
    template<typename TypeTag>
    class EbosProblemFlow: public EbosProblem<TypeTag>{
    public:
        using Parent = EbosProblem<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using TimeStepper =  AdaptiveTimeSteppingEbos<TypeTag>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        EbosProblemFlow(Simulator& simulator): EbosProblem<TypeTag>(simulator){
        }
        void timeIntegration()
        {
            if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Start TimeIntegration-------------------\n"
                << std::flush;                                                     
            }
            Parent::timeIntegration();
        }
        void finishInit(){
            Parent::finishInit();
            this->simulator().model().addAuxiliaryModule(&monitorAux_);
        }
    private:
        using MonitorAuxType = MonitoringAuxModule<TypeTag>;
        MonitorAuxType monitorAux_;    
    
        //private:
        //std::unique_ptr<TimeStepper> adaptiveTimeStepping_;
    };

    template<typename TypeTag>
    class BlackoilWellModelFvExtra: public BlackoilWellModel<TypeTag>{
        using Parent = BlackoilWellModel<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    public:
        BlackoilWellModelFvExtra(Simulator& ebosSimulator): Parent(ebosSimulator) {};
        void beginIteration(){
            Parent::beginIteration();
            std::cout << "EclWellModelFvExtra begin iteration" << std::endl;
        }
        void endIteration(){
            Parent::endIteration();
            std::cout << "EclWellModelFvExtra end iteration" << std::endl;
        }
    };
        

    
    template<typename TypeTag>
    class EclNewtonMethodLinesearch: public EclNewtonMethod<TypeTag>{
    public:
        using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
        using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
        using Parent = EclNewtonMethod<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using TimeStepper =  AdaptiveTimeSteppingEbos<TypeTag>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        EclNewtonMethodLinesearch(Simulator& simulator): Parent(simulator){
        }

        void update_(SolutionVector& nextSolution,
                     const SolutionVector& currentSolution,
                     const GlobalEqVector& solutionUpdate,
                     const GlobalEqVector& currentResidual){
            //if (this->gridView().comm().rank() == 0){
                std::cout << "----------------------Newton update_ With Linesearch------------\n"
                          << std::flush;
                //}
            auto& linearizer = this->model().linearizer();
            auto& residual = linearizer.residual();//use reference
            auto newSolutionUpdate = solutionUpdate;
            Scalar error = this->norm(residual);
            Scalar well_error = this->well_error();
            int iter = 0;
            std::cout << "Linesarch itr " << iter
                      << " ResError " << error
                      << " WellError " << well_error << std::endl;    
            Parent::update_(nextSolution,currentSolution,solutionUpdate,currentResidual);
            this->problem().endIteration();               
            Scalar target = error;
            error = 1e99;
            while((error >= target) && (iter < 5)){
                iter +=1;
                // update the solution in aux modules
                //this->beginIteration_();
                this->problem().beginIteration();
                this->linearizeDomain_();
                this->linearizeAuxiliaryEquations_();
                this->preSolve_(nextSolution,residual);
                this->postSolve_(currentSolution,
                                 currentResidual,
                                 newSolutionUpdate);
                //this->endIteration_(nextSolution, currentSolution);                
                error = this->norm(residual);
                well_error = this->well_error();
                std::cout << "Linesarch itr " << iter
                          << " ResError " << error
                          << " WellError " << well_error << std::endl;
                newSolutionUpdate *= 0.5;
                Parent::update_(nextSolution,currentSolution,newSolutionUpdate,currentResidual);
                this->problem().endIteration();               
            }
            //previousError_ = error;
        }
        Scalar well_error() const{
            Scalar sum = 0;
            const auto& well_container = this->problem().wellModel().localNonshutWells();
            using StdWell = StandardWell<TypeTag>;
            for (const auto& well: well_container) {
                StdWell* stdwell = dynamic_cast<StdWell*>(well.get());
                Scalar well_norm;
                if(!stdwell==0){
                    auto& residual = stdwell->linSys().residual();
                    well_norm = this->norm(residual);
                }else{
                    throw std::runtime_error("not standard wells");
                }                                    
                sum += well_norm;                        
            }
            return sum;
        }

        template<class VectorVector>
        Scalar norm(const VectorVector& vec) const{
            Scalar sum = 0;
            for(const auto& bvec : vec){
                for(const auto& val: bvec){
                    sum += val*val;
                }
            }
            // maybe add rhs of auxModules i.e. wells and aquifers
            return sum;
        }
        
        bool proceed_(){
            std::cout << "----------------------Newton with early exit-------------------\n"
                      << std::flush;                
            bool proceed = Parent::proceed_();
            return proceed;
        }
        //private:
        //Scalar previousError_;
        
    };
    
} // namespace Opm

// the current code use eclnewtonmethod adding other conditions to proceed_ should do the trick for KA
// adding linearshe sould be chaning the update_ function in the same class with condition that the error is reduced.
// the trick is to be able to recalculate the residual from here.
// unsure where the timestepping is done from suggestedtime??
// suggestTimeStep is taken from newton solver in problem.limitTimestep
namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowProblemEbos {
    using InheritsFrom = std::tuple<EbosTypeTag>;
};
}

// Set the problem class
template<class TypeTag>
struct Problem<TypeTag, TTag::EclFlowProblemEbos> {
    using type = EbosProblemFlow<TypeTag>;
};

template<class TypeTag>
struct Model<TypeTag, TTag::EclFlowProblemEbos> {
    using type = BlackOilModelFvNoCache<TypeTag>;
    // using type = BlackOilModelFvLocal<TypeTag>;
};    

template<class TypeTag>
struct EclWellModel<TypeTag, TTag::EclFlowProblemEbos> {
    using type = BlackoilWellModelFvExtra<TypeTag>;
};
    
template<class TypeTag>
struct NewtonMethod<TypeTag, TTag::EclFlowProblemEbos> {
    using type = EclNewtonMethodLinesearch<TypeTag>;
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

template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemEbos> {
    //using type = EclBlackOilIntensiveQuantities<TypeTag>;
    using type = BlackOilIntensiveQuantitiesSimple<TypeTag>;
    //using type = BlackOilIntensiveQuantities<TypeTag>;
    //using type = BlackOilIntensiveQuantitiesDryGas<TypeTag>;
};

// template<class TypeTag>
// struct Linearizer<TypeTag, TTag::EclFlowProblemEbos> { using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::EclFlowProblemEbos> { using type = BlackOilLocalResidualTPFA<TypeTag>; };
    
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::EclFlowProblemEbos> { static constexpr bool value = false; };

template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::EclFlowProblemEbos> { static constexpr bool value = false; };    

//static constexpr bool has_disgas_in_water = getPropValue<TypeTag, Properties::EnableDisgasInWater>();    
    
template<class TypeTag>
struct Simulator<TypeTag, TTag::EclFlowProblemEbos> { using type = Opm::Simulator<TypeTag>; };    
// // Set the problem class
// template<class TypeTag>
// struct Problem<TypeTag, TTag::EbosTypeTag> {
//     using type = EbosProblem<TypeTag>;
// };
    
template<class TypeTag>
struct EclEnableAquifers<TypeTag, TTag::EclFlowProblemEbos> {
     static constexpr bool value = false;
};
}
}
int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::EclFlowProblemEbos;
    Opm::registerEclTimeSteppingParameters<TypeTag>();
    return Opm::start<TypeTag>(argc, argv);
}
