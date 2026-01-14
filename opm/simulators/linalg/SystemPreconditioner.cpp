#include <config.h>
#include <opm/simulators/linalg/SystemPreconditioner.hpp>
#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>



namespace Opm {
    SystemPreconditioner::SystemPreconditioner(const SystemMatrix& S,const std::function<ResVector()> &weightsCalculator,int pressureIndex, 
    const Opm::PropertyTree& prm) :
        S_(S), prm_(prm) {
            auto rop = std::make_unique<ResOperator>(S[_0][_0]);
            auto wop = std::make_unique<WellOperator>(S[_1][_1]);
            auto resprm = prm_.get_child("reservoir_solver");
            auto wellprm = prm_.get_child("well_solver");
            //std::function<ResVector()> weightsCalculatorRes;
            auto rsol = std::make_unique<ResFlexibleSolverType>(*rop, resprm,
                                                             weightsCalculator,
                                                             pressureIndex);
            std::function<WellVector()> weightsCalculatorWell;
            auto wsol = std::make_unique<WellFlexibleSolverType>(*wop, wellprm,
                                                             weightsCalculatorWell,
                                                             pressureIndex); 
            this->rop_ = std::move(rop);
            this->wop_ = std::move(wop);
            this->resSolver_ = std::move(rsol);
            this->wellSolver_ = std::move(wsol);

        }
    
    void 
    SystemPreconditioner::apply(SystemVector& v, const SystemVector& d) {
            // change order?
            Dune::InverseOperatorResult well_result;
            const auto& C = S_[_1][_0];
            const auto& B = S_[_0][_1];
            const auto& r_r = d[_0];
            const auto& r_w = d[_1];

            auto wRes = r_w;
            C.mv(v[_0], wRes);
            WellVector wSol(wRes.size());
            //wSol = v[_1];
            wSol = 0.0;
            double well_tol = prm_.get<double>("well_solver.tol");
            wellSolver_->apply(wSol, wRes, well_tol, well_result);
            v[_1] = wSol;
            // Use reservoir solver
            
            auto resRes = r_r;
            ResVector resSol(resRes.size());
            resSol = 0.0;
            B.mmv(wSol, resRes);
            Dune::InverseOperatorResult res_result;
            double res_tol = prm_.get<double>("reservoir_solver.tol");
            resSolver_->apply(resSol, resRes, res_tol, res_result);
            v[_0] = resSol;
            // solve well again
            
            wRes = d[_1];
            C.mmv(v[_0], wRes);
            wSol = 0.0;
            wellSolver_->apply(wSol, wRes, well_tol, well_result);
            v[_1] = wSol;
            
            // update well rhs

            // update residual
            // S_[1][_0].mv(vPart, wellPart, -1.0, 1.0);
        
    }

}
const int numResDofs = 3;
const int numWellDofs = 4;
using CommSeq = Dune::Amg::SequentialInformation;
using CommPar = Dune::OwnerOverlapCopyCommunication<int, int>;
using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
using WellVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
using WellOp = Dune::MatrixAdapter<WWMatrix, WellVector, WellVector>;
using WellFlexibleSolverType = Dune::FlexibleSolver<WellOp>;  
template class Dune::FlexibleSolver<WellOp>;
//template class Dune::FlexibleSolver<WellOp>;
/* template Dune::FlexibleSolver<WellOp>::FlexibleSolver(WellOp& op,
                     const CommPar& comm,
                     const Opm::PropertyTree& prm,
                     const std::function<WellVector()>& weightsCalculator,
                     std::size_t pressureIndex); */
//template class Opm::PreconditionerFactory<WellOp, CommSeq>;
//template class Opm::PreconditionerFactory<WellOp, CommPar>;                                          
// Define matrix and vector typ