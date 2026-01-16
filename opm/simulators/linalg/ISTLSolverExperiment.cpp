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
#include <opm/common/ErrorMacros.hpp>
namespace Dune
{
    class SeqComm
    {
    public:
        using size_type = std::size_t;
        // SolverCategory::Category category () const {
        //     return category_;
        // }

        // const Communication<MPI_Comm>& communicator() const
        // {
        //   return cc;
        // }
        static constexpr size_type size()
        {
            return 0;
        }
        template <typename T>
        void project(T & /*x*/) const
        {
            // No operation for sequential communicator
        }
        template <typename T1, typename T2>
        void dot(const T1 &x, const T1 &y, T2 &result) const
        {
            result = x.dot(y);
        }
        template <typename T>
        double norm(const T &x) const
        {
            return x.two_norm();
        };
        template <typename T>
        void copyOwnerToAll(const T &x,T &y) const
        {
            y = x;
        }
    };

  //template <typename ColCom> Communication<MPI_Comm>
    class JacComm
    {
    public:
        using size_type = std::size_t;
        // SolverCategory::Category category () const {
        //     return category_;
        // }

        // const Communication<MPI_Comm>& communicator() const
        // {
        //   return cc;
        // }
    public:
      //JacComm(const ColCom& colcom): colcom_(colcom){}
      JacComm(): colcom_(MPI_COMM_WORLD){}
        static constexpr size_type size()
        {
            return 0;
        }
        template <typename T>
        void project(T & /*x*/) const
        {
            // No operation for sequential communicator
        }
      
        template <typename T1, typename T2>
        void dot(const T1 &x, const T1 &y, T2 &result) const
        {
            result = x.dot(y);
            result = colcom_.sum(result);
        }
      
        template <typename T>
        double norm(const T &x) const
        {
            double result = x.two_norm();
            result = colcom_.sum(result);
            return result;
        }
      
        template <typename T>
        void copyOwnerToAll(const T &x,T &y) const
        {
            y = x;
        }
       private:
       Communication<MPI_Comm> colcom_;
    };
  
    template <typename... Args>
    class MultiCommunicator
        : public std::tuple<Args...>
    {
        /** \brief Helper type */
        typedef std::tuple<Args...> TupleType;
        typedef MultiCommunicator<Args...> type;
        using field_type = double;

    public:
        using std::tuple<Args...>::tuple;
        using size_type = std::size_t;

            /**
     * @brief Get Solver Category.
     * @return The Solver Category.
     */
        // SolverCategory::Category category () const {
        //     return category_;
        // }

        // const Communication<MPI_Comm>& communicator() const
        // {
        //   return cc;
        // }

        static constexpr size_type size()
        {
            return sizeof...(Args);
        }
        /** \brief Number of elements
         */
        static constexpr size_type N()
        {
            return sizeof...(Args);
        }
        template <size_type index>
        typename std::tuple_element<index, TupleType>::type &
        operator[]([[maybe_unused]] const std::integral_constant<size_type, index> indexVariable)
        {
            return std::get<index>(*this);
        }
        template <size_type index>
        const typename std::tuple_element<index, TupleType>::type &
        operator[]([[maybe_unused]] const std::integral_constant<size_type, index> indexVariable) const
        {
            return std::get<index>(*this);
        }

        template <typename T>
        void project(T &x) const
        {
            using namespace Dune::Hybrid;
            //auto size = index_constant<sizeof...(Args)>();
            forEach(integralRange(Hybrid::size(*this)), [&](auto &&i)
                    { (*this)[i].project(x[i]); });
        }
        template <typename T1, typename T2>
        void dot(const T1 &x, const T1 &y, T2& result) const
        {
            result = field_type(0);
            using namespace Dune::Hybrid;
            forEach(integralRange(Hybrid::size(*this)), [&](auto &&i)
            { double result_tmp = 0;
              (*this)[i].dot(x[i],y[i],result_tmp);
              result += result_tmp;
              //std::cout << " Dot partial result " << i << " r " << result <<std::endl;
            } );
        }
        template <typename T>
        field_type norm(const T &x) const
        {
            using namespace Dune::Hybrid;
            return accumulate(integralRange(Hybrid::size(*this)), field_type(0), [&](auto &&a, auto &&i)
                              { return a + (*this)[i].norm(x[i]); });
        }
        template <typename T>
        void copyOwnerToAll(const T &x,T &y) const
        {
            using namespace Dune::Hybrid;
            forEach(integralRange(Hybrid::size(*this)), [&](auto &&i)
                    { (*this)[i].copyOwnerToAll(x[i], y[i]); });
        }
        // Explicit template instantiations
    };
}
namespace Opm
{
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

    Dune::InverseOperatorResult solveSystem(const SystemMatrix& S, SystemVector& x, const SystemVector& b, 
        const std::function<RVector()> &weightCalculator,int pressureIndex, const Opm::PropertyTree& prm)
    {
        // Here we would implement the solver logic for the system S * x = b
        // This is a placeholder implementation
        int verbosity = prm.get<int>("verbosity");           // Reduce output verbosity
        if(verbosity){
            std::cout << "Solving system with merged matrices..." << std::endl;
        }
         const Dune::MatrixAdapter<SystemMatrix, SystemVector, SystemVector> S_linop(S);
         //const Dune::OverlappingSchwarzOperator<SystemVector, SystemVector, Dune::OwnerOverlapCopyCommunication<int, int> > S_linop(S_linop, comm);
    //TailoredPrecondDiag precond(S,prm);
    Opm::PropertyTree precond_prm = prm.get_child("preconditioner");
    SystemPreconditioner precond(S,weightCalculator, pressureIndex, precond_prm);
    
    // Set solver parameters
    double linsolve_tol = prm.get<double>("tol");  // Less strict tolerance
    int max_iter = prm.get<int>("maxiter");           // Limit iterations
   
    Dune::InverseOperatorResult result;
    // Create and run the solver with error handling
    try {
        if(verbosity > 0){
        std::cout << "Solving system with BiCGSTAB solver..." << std::endl;
        }
        const std::string solver_type = prm.get<std::string>("solver");
        using AbstractSolverType = Dune::InverseOperator<SystemVector, SystemVector>;
        std::shared_ptr<AbstractSolverType> linsolver;
        if( solver_type == "bicgstab"){
            linsolver = std::make_shared<Dune::BiCGSTABSolver<SystemVector>>(
                                                                             S_linop,
                                                                             precond,
                                                                             linsolve_tol,
                                                                             max_iter,
                                                                             verbosity
                                                                             );
        } else if ( solver_type == "fgmres"){
          int restart = prm.get<int>("restart", 15);
          linsolver = std::make_shared<Dune::RestartedGMResSolver<SystemVector>>(
                                                                               S_linop,
                                                                               precond,
                                                                               linsolve_tol,
                                                                               restart,
                                                                               max_iter,
                                                                               verbosity
                                                                             );
        }else {
          OPM_THROW(std::invalid_argument,
                      "Properties: Solver " + solver_type + " not known.");
        }
        auto residual(b);
        linsolver->apply(x, residual, result);
        //assert(false);//debug
        // Print results
        if(verbosity > 10){
        std::cout << "\nSolver results:" << std::endl;
        std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << std::endl;
        std::cout << "  Iterations: " << result.iterations << std::endl;
        std::cout << "  Reduction: " << result.reduction << std::endl;
        std::cout << "  Elapsed time: " << result.elapsed << " seconds" << std::endl;
        std::cout << "\nMatrix merger example completed successfully!" << std::endl;
        }
        //printVector("Solution x_r", x[_0]);
        //printVector("Solution x_w", x[_1]);
    }
    catch (const std::exception& e) {
        std::cerr << "Error solving system: " << e.what() << std::endl;
        //return result;
    }
    
    
    return result;
    }

    Dune::InverseOperatorResult solveSystem(const SystemMatrix& S, SystemVector& x, const SystemVector& b, 
        const std::function<RVector()> &weightCalculator,int pressureIndex, const Opm::PropertyTree& prm,
        const Dune::OwnerOverlapCopyCommunication<int, int>& comm)
    {
        // Here we would implement the solver logic for the system S * x = b
        // This is a placeholder implementation
          const bool is_iorank = comm.communicator().rank() == 0;
        const int verbosity = is_iorank ? prm.get<int>("verbosity", 0) : 0;
        if(verbosity){
            std::cout << "Solving system with merged matrices..." << std::endl;
        }
         //const Dune::MatrixAdapter<SystemMatrix, SystemVector, SystemVector> S_linop(S);
        using WellComm = Dune::JacComm;
        using SystemComm = Dune::MultiCommunicator<const Dune::OwnerOverlapCopyCommunication<int, int>&,const WellComm&>;
         WellComm seqComm;
         SystemComm systemComm(comm,seqComm);
         
         const Dune::OverlappingSchwarzOperator<SystemMatrix, SystemVector, SystemVector, SystemComm > S_linop(S, systemComm);
         std::shared_ptr< Dune::ScalarProduct<SystemVector> > scalarproduct = Dune::createScalarProduct<SystemVector, SystemComm>(systemComm, S_linop.category());
         //
    //TailoredPrecondDiag precond(S,prm);
    Opm::PropertyTree precond_prm = prm.get_child("preconditioner");
    SystemPreconditioner precond(S,weightCalculator, pressureIndex, precond_prm);
    Dune::BlockPreconditioner<SystemVector, SystemVector, SystemComm,SystemPreconditioner> sysprecond(precond, systemComm);
    
    // Set solver parameters
    double linsolve_tol = prm.get<double>("tol");  // Less strict tolerance
    int max_iter = prm.get<int>("maxiter");           // Limit iterations
   
    Dune::InverseOperatorResult result;
    // Create and run the solver with error handling
    try {
        if(verbosity > 0){
        std::cout << "Solving system with BiCGSTAB solver parallel...rank.." << comm.communicator().rank() << std::endl;
        }
        using AbstractSolverType = Dune::InverseOperator<SystemVector, SystemVector>;
        std::shared_ptr<AbstractSolverType> linsolver;
        const std::string solver_type = prm.get<std::string>("solver");
        if( solver_type == "bicgstab"){
            linsolver = std::make_shared<Dune::BiCGSTABSolver<SystemVector>>(
                                                                             S_linop,
                                                                             *scalarproduct,
                                                                             sysprecond,
                                                                             linsolve_tol,
                                                                             max_iter,
                                                                             verbosity
                                                                             );
        }else if ( solver_type == "fgmres"){
          int restart = prm.get<int>("restart", 15);
          linsolver = std::make_shared<Dune::RestartedGMResSolver<SystemVector>>(
                                                                             S_linop,
                                                                             *scalarproduct,
                                                                             sysprecond,
                                                                             linsolve_tol,
                                                                             restart,
                                                                             max_iter,
                                                                             verbosity
                                                                             );
        }else {
          OPM_THROW(std::invalid_argument,
                      "Properties: Solver " + solver_type + " not known.");
        }
        auto residual(b);
        linsolver->apply(x, residual, result);
        //assert(false);
                      
        // Print results
        if(verbosity > 10){
        std::cout << "\nSolver results:" << std::endl;
        std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << std::endl;
        std::cout << "  Iterations: " << result.iterations << std::endl;
        std::cout << "  Reduction: " << result.reduction << std::endl;
        std::cout << "  Elapsed time: " << result.elapsed << " seconds" << std::endl;
        std::cout << "\nMatrix merger example completed successfully!" << std::endl;
        }
        //printVector("Solution x_r", x[_0]);
        //printVector("Solution x_w", x[_1]);
    }
    catch (const std::exception& e) {
        std::cerr << "Error solving system: " << e.what() << std::endl;
        //return result;
    }
    
    
    return result;
    }
}
    }
