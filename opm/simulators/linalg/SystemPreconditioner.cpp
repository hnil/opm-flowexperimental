#include <config.h>
#include <opm/simulators/linalg/SystemPreconditioner.hpp>
#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>

namespace Opm {
    SystemPreconditioner::SystemPreconditioner(const SystemMatrix& S, Opm::PropertyTree& prm) :
        A_diag_(diagvec(S[_0][_0])), M_diag_(diagvec(S[_1][_1])), prm_(prm) {
            auto rop = std::make_unique<ResOperator>(S[_0][_0]);
            auto wop = std::make_unique<WellOperator>(S[_1][_1]);
            auto resprm = prm_.get_child("reservoir_solver");
            auto wellprm = prm_.get_child("well_solver");
            std::function<ResVector()> weightsCalculatorRes;
            int pressureIndex = 0; // Assuming pressure is the first variable
            auto rsol = std::make_unique<ResFlexibleSolverType>(*rop, resprm,
                                                             weightsCalculatorRes,
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
        bool simple = false;
        v = d;
        //return;
        if (simple)
        {
            for (size_t i = 0; i != A_diag_.size(); ++i)
            {
                for (size_t j = 0; j < numResDofs; ++j)
                {
                    v[_0][i][j] = d[_0][i][j] / A_diag_[i][j];
                }
            }
            for (size_t i = 0; i != M_diag_.size(); ++i)
            {
                for (size_t j = 0; j < numWellDofs; ++j)
                {
                    //v[_1][i][j] = d[_1][i][j] / M_diag_[i][j];
                }
            }
        }
        else
        {
            // change order?
            Dune::InverseOperatorResult well_result;
            WellVector wellPart = d[_1];
            WellVector vWellPart(wellPart.size());
            vWellPart = wellPart;
            double well_tol = prm_.get<double>("well_solver.tol");
            wellSolver_->apply(vWellPart, wellPart, well_tol, well_result);
            v[_1] = vWellPart;
            // Use reservoir solver
            ResVector resPart = d[_0];
            // updete residual
            // S_[0][_1].mv(vWellPart, resPart, -1.0, 1.0);
            Dune::InverseOperatorResult res_result;
            ResVector vPart(resPart.size());
            double res_tol = prm_.get<double>("reservoir_solver.tol");
            resSolver_->apply(vPart, resPart, 1e-3, res_result);
            v[_0] = vPart;
            // update well rhs

            // update residual
            // S_[1][_0].mv(vPart, wellPart, -1.0, 1.0);
        }
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