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
#include <ebos/eclproblem.hh>
#include <ebos/eclnewtonmethod.hh>
#include <ebos/ebos.hh>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/models/blackoil/blackoilintensivequantitiessimple.hh>
namespace Opm{    
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
        
        //private:
        //std::unique_ptr<TimeStepper> adaptiveTimeStepping_;
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
        const IntensiveQuantities* cachedIntensiveQuantities(unsigned globalIdx, unsigned timeIdx) const{
            //std::cout << "----------------------Update cached quantites with globalIdx-------------------\n";
            //const IntensiveQuantities* intQuants = Parent::cachedIntensiveQuantities(globalIdx,timeIdx);return intQuants;
            //const auto& primaryVar = this->solution(timeIdx)[globalIdx];
            //const auto& problem = this->simulator_.problem();
            IntensiveQuantities* intQuants = &(this->intensiveQuantityCache_[timeIdx][globalIdx]);
            if (!(this->enableIntensiveQuantityCache_) ||
                !(this->intensiveQuantityCacheUpToDate_[timeIdx][globalIdx])){
                ///intQuants->update(problem,primaryVar, globalIdx, timeIdx);
                return 0;
            }else{
                return intQuants;
            }
            
        }
        IntensiveQuantities intensiveQuantities(unsigned globalIdx, unsigned timeIdx) const{
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
                
            Parent::update_(nextSolution,currentSolution,solutionUpdate,currentResidual);
            // calculate new residual
            Scalar previousError = this->norm(currentResidual);
            Scalar well_errors = this->well_error();
            int iter = 0;
            std::cout << "Previous error" << previousError << std::endl;
            Scalar target = previousError;
            auto& linearizer = this->model().linearizer();
            auto& residual = linearizer.residual();//use reference
            auto newSolutionUpdate = solutionUpdate;
            Scalar error = this->norm(residual);
            std::cout << "Linesarch itr org" << iter << " error " << error << std::endl;
            while((error >= target) && (iter < 5)){
                iter +=1;
                newSolutionUpdate *= 0.0;
                // update the solution in aux modules
                
                
                //this->beginIteration_();
                //problem().beginIteration();
                residual = 0;
                this->linearizeDomain_();
                //this->linearizeAuxiliaryEquations_();
                //this->preSolve_(nextSolution,residual);
                error = norm<Scalar>(residual);
                //this->postSolve_(currentSolution,
                //                 currentResidual,
                //                 newSolutionUpdate);
                //Parent::update_(nextSolution,currentSolution,newSolutionUpdate,currentResidual);
                
                //this->endIteration_(nextSolution, currentSolution);
                //problem().endIterations();
                
                std::cout << "Linesarch itr " << iter << " error " << error << std::endl;
            }
            //previousError_ = error;
        }
        Scalar well_error() const{
            Scalar sum = 0;
            const auto& well_container = this->localNonShutWells(); 
            for (const auto& well: well_container) {
                StandardWell* stdwell = dynamic_cast<StandardWell>(well);
                Scalar well_norm;
                if(!stdwell==0){
                    auto& residual = stdwell->residual();
                    well_norm = this->norm(residual);
                }else{
                    std::error("not standard wells");
                }                                    
                sum += well_norm;                        
            }
            return sum;
        }

        
        Scalar norm(const GlobalEqVector& vec) const{
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
    using type = BlackOilModelFv<TypeTag>;
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
    
// template<class TypeTag>
// struct EclEnableAquifers<TypeTag, TTag::EclFlowProblemTest> {
//     static constexpr bool value = false;
// };
}
}
int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::EclFlowProblemEbos;
    return Opm::start<TypeTag>(argc, argv);
}
