#include "config.h"
#include <iostream>
#include <dune/istl/bvector.hh>
#include <opm/simulators/linalg/matrixblock.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include "WellMatrixMerger.hpp"
#include "SystemPreconditioner.hpp"
namespace Opm {
        namespace SystemSolver {
  const int numResDofs = 3;
            const int numWellDofs = 4;

                // Define matrix and vector types
             //using RRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numResDofs>>;
             using RRMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, numResDofs, numResDofs>>;
                //using RWtype = Dune::FieldMatrix<double, numResDofs, numWellDofs>;
                using RWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numWellDofs>>;
                using WRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numResDofs>>;
                using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
                using RVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
                using WVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
                  using SystemMatrix = Dune::MultiTypeBlockMatrix<
                Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
                Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>
        >;

        using SystemMatrix = Dune::MultiTypeBlockMatrix<
                Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
                Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>
        >; 
        using SystemVector = Dune::MultiTypeBlockVector<RVector, WVector>;

    void solveSystem(const SystemMatrix& S, SystemVector& x, const SystemVector& b,Opm::PropertyTree& prm)
    {
        // Here we would implement the solver logic for the system S * x = b
        // This is a placeholder implementation
        std::cout << "Solving system with merged matrices..." << std::endl;
         const Dune::MatrixAdapter<SystemMatrix, SystemVector, SystemVector> S_linop(S);
    //TailoredPrecondDiag precond(S,prm);
    Opm::PropertyTree precond_prm = prm.get_child("preconditioner");
    SystemPreconditioner precond(S,precond_prm);
    
    // Set solver parameters
    double linsolve_tol = prm.get<double>("tol");  // Less strict tolerance
    int max_iter = prm.get<int>("maxiter");           // Limit iterations
    int verbosity = prm.get<int>("verbosity");           // Reduce output verbosity
    
    // Create and run the solver with error handling
    try {
        std::cout << "Solving system with BiCGSTAB solver..." << std::endl;
        
        auto solver = Dune::BiCGSTABSolver<SystemVector>(
            S_linop,
            precond,
            linsolve_tol,
            max_iter,
            verbosity
        );
        
        Dune::InverseOperatorResult result;
        auto residual(b);
        solver.apply(x, residual, result);
        
        // Print results
        std::cout << "\nSolver results:" << std::endl;
        std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << std::endl;
        std::cout << "  Iterations: " << result.iterations << std::endl;
        std::cout << "  Reduction: " << result.reduction << std::endl;
        std::cout << "  Elapsed time: " << result.elapsed << " seconds" << std::endl;
        
        //printVector("Solution x_r", x[_0]);
        //printVector("Solution x_w", x[_1]);
    }
    catch (const std::exception& e) {
        std::cerr << "Error solving system: " << e.what() << std::endl;
    }
    
    std::cout << "\nMatrix merger example completed successfully!" << std::endl;
    
    }
}
}

