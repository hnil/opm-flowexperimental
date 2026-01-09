/*
  Copyright 2023 OPM Contributors

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

#include <config.h>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <map>
#include <string>
#include <functional>

// Define constants
const int numResDofs = 3;
const int numWellDofs = 4;

// Define matrix and vector types
using RRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numResDofs>>;
using RWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numWellDofs>>;
using WRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numResDofs>>;
using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
using RVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
using WVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;

// Define system matrix and vector types
using SystemMatrix = Dune::MultiTypeBlockMatrix<
    Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
    Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>
>;
using SystemVector = Dune::MultiTypeBlockVector<RVector, WVector>;

// Define helper constants for accessing MultiTypeBlockMatrix elements
constexpr auto _0 = Dune::Indices::_0;
constexpr auto _1 = Dune::Indices::_1;

// Helper function to create a matrix with a specific pattern and values
template<class Matrix, class Block>
void createMatrix(Matrix& matrix, int numRows, int numCols, const std::vector<int>& colIndices, const Block& blockValue)
{
    matrix.setSize(numRows, numCols);
    matrix.setBuildMode(Matrix::row_wise);
    
    // Define the sparsity pattern
    for (auto row = matrix.createbegin(); row != matrix.createend(); ++row) {
        for (const auto& col : colIndices) {
            if (col < numCols) {
                row.insert(col);
            }
        }
    }
    
    // Set the values
    for (int row = 0; row < numRows; ++row) {
        for (const auto& col : colIndices) {
            if (col < numCols) {
                matrix[row][col] = blockValue;
            }
        }
    }
}

// Helper function to print a matrix
template<class Matrix>
void printMatrix(const std::string& name, const Matrix& matrix)
{
    std::cout << "Matrix " << name << " (" << matrix.N() << " x " << matrix.M() << "):" << std::endl;
    
    for (auto row = matrix.begin(); row != matrix.end(); ++row) {
        std::cout << "  Row " << row.index() << ": ";
        for (auto col = row->begin(); col != row->end(); ++col) {
            std::cout << "Col " << col.index() << " ";
            // Print the first element of each block
            std::cout << "(" << (*col)[0][0] << ") ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Helper function to print a vector
template<class Vector>
void printVector(const std::string& name, const Vector& vector)
{
    std::cout << "Vector " << name << " (size " << vector.size() << "):" << std::endl;
    std::cout << "  [";
    for (size_t i = 0; i < vector.size(); ++i) {
        std::cout << "(";
        for (size_t j = 0; j < vector[i].size(); ++j) {
            std::cout << std::fixed << std::setprecision(4) << vector[i][j];
            if (j < vector[i].size() - 1) std::cout << ", ";
        }
        std::cout << ")";
        if (i < vector.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl << std::endl;
}

// Extract diagonal elements from a matrix
template<class Matrix>
auto diagvec(const Matrix& matrix) {
    using BlockType = typename Matrix::block_type;
    using FieldType = typename BlockType::field_type;
    using VectorBlock = Dune::FieldVector<FieldType, BlockType::rows>;
    using ResVector = Dune::BlockVector<VectorBlock>;
    
    ResVector diag(matrix.N());
    for (size_t i = 0; i < matrix.N(); ++i) {
        auto row = matrix[i];
        auto col = row.find(i);
        if (col != row.end()) {
            // Extract diagonal elements
            for (size_t j = 0; j < BlockType::rows; ++j) {
                diag[i][j] = (*col)[j][j];
            }
        } else {
            // If diagonal element doesn't exist, set to 1.0
            for (size_t j = 0; j < BlockType::rows; ++j) {
                diag[i][j] = 1.0;
            }
        }
    }
    return diag;
}

// Simple diagonal preconditioner for the system
class TailoredPrecondDiag : public Dune::Preconditioner<SystemVector, SystemVector>
{
public:
    using ResVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
    using WellVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
    
    TailoredPrecondDiag(const SystemMatrix& S) :
        A_diag_(diagvec(S[_0][_0])), M_diag_(diagvec(S[_1][_1])) {}
    
    virtual void apply(SystemVector& v, const SystemVector& d) override {
        for (size_t i = 0; i != A_diag_.size(); ++i) {
            for (size_t j = 0; j < numResDofs; ++j) {
                v[_0][i][j] = d[_0][i][j] / A_diag_[i][j];
            }
        }
        for (size_t i = 0; i != M_diag_.size(); ++i) {
            for (size_t j = 0; j < numWellDofs; ++j) {
                v[_1][i][j] = d[_1][i][j] / M_diag_[i][j];
            }
        }
    }
    
    virtual void pre(SystemVector& x, SystemVector& b) override {}
    virtual void post(SystemVector& x) override {}
    virtual Dune::SolverCategory::Category category() const override {
        return Dune::SolverCategory::sequential;
    }
    
private:
    const ResVector A_diag_;
    const WellVector M_diag_;
};

// Count non-zeros in a matrix row
template<class Row>
int countRowEntries(const Row& row) {
    return std::distance(row.begin(), row.end());
}

// Merge B matrices (well-to-reservoir) by vertical stacking
void mergeWRMatrices(const std::vector<WRMatrix>& matrices, WRMatrix& mergedMatrix)
{
    if (matrices.empty()) {
        return;
    }
    
    // Count total rows and get number of columns
    int numCols = matrices[0].M();
    int totalRows = 0;
    for (const auto& matrix : matrices) {
        totalRows += matrix.N();
    }
    
    // Set up the merged matrix with random build mode (required for setrowsize)
    mergedMatrix.setSize(totalRows, numCols);
    mergedMatrix.setBuildMode(WRMatrix::random);
    
    // Phase 1: Set row sizes
    int rowOffset = 0;
    for (const auto& matrix : matrices) {
        for (int i = 0; i < static_cast<int>(matrix.N()); ++i) {
            int nnz = countRowEntries(matrix[i]);
            mergedMatrix.setrowsize(rowOffset + i, nnz);
        }
        rowOffset += matrix.N();
    }
    mergedMatrix.endrowsizes();
    
    // Phase 2: Set column indices
    rowOffset = 0;
    for (const auto& matrix : matrices) {
        for (int i = 0; i < static_cast<int>(matrix.N()); ++i) {
            for (auto colIt = matrix[i].begin(); colIt != matrix[i].end(); ++colIt) {
                mergedMatrix.addindex(rowOffset + i, colIt.index());
            }
        }
        rowOffset += matrix.N();
    }
    mergedMatrix.endindices();
    
    // Phase 3: Copy values
    rowOffset = 0;
    for (const auto& matrix : matrices) {
        for (int i = 0; i < static_cast<int>(matrix.N()); ++i) {
            for (auto colIt = matrix[i].begin(); colIt != matrix[i].end(); ++colIt) {
                mergedMatrix[rowOffset + i][colIt.index()] = *colIt;
            }
        }
        rowOffset += matrix.N();
    }
}

// Merge C matrices (reservoir-to-well) by horizontal stacking
void mergeRWMatrices(const std::vector<RWMatrix>& matrices, RWMatrix& mergedMatrix)
{
    if (matrices.empty()) {
        return;
    }
    
    // All matrices should have the same number of rows (reservoir cells)
    int numRows = matrices[0].N();
    
    // Count total columns (wells)
    int totalCols = 0;
    for (const auto& matrix : matrices) {
        totalCols += matrix.M();
    }
    
    // Set up the merged matrix with random build mode (required for setrowsize)
    mergedMatrix.setSize(numRows, totalCols);
    mergedMatrix.setBuildMode(RWMatrix::random);
    
    // Phase 1: Set row sizes
    for (int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
        int nnz = 0;
        for (const auto& matrix : matrices) {
            if (rowIdx < static_cast<int>(matrix.N())) {
                nnz += countRowEntries(matrix[rowIdx]);
            }
        }
        mergedMatrix.setrowsize(rowIdx, nnz);
    }
    mergedMatrix.endrowsizes();
    
    // Phase 2: Set column indices
    for (int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
        int colOffset = 0;
        for (const auto& matrix : matrices) {
            if (rowIdx < static_cast<int>(matrix.N())) {
                for (auto colIt = matrix[rowIdx].begin(); colIt != matrix[rowIdx].end(); ++colIt) {
                    mergedMatrix.addindex(rowIdx, colOffset + colIt.index());
                }
            }
            colOffset += matrix.M();
        }
    }
    mergedMatrix.endindices();
    
    // Phase 3: Copy values
    int colOffset = 0;
    for (const auto& matrix : matrices) {
        for (int i = 0; i < static_cast<int>(matrix.N()); ++i) {
            for (auto colIt = matrix[i].begin(); colIt != matrix[i].end(); ++colIt) {
                mergedMatrix[i][colOffset + colIt.index()] = *colIt;
            }
        }
        colOffset += matrix.M();
    }
}

// Create diagonal block matrix from D matrices (well-to-well)
void createDiagonalBlockMatrix(const std::vector<WWMatrix>& matrices, WWMatrix& mergedMatrix)
{
    if (matrices.empty()) {
        return;
    }
    
    // Count total size
    int totalSize = 0;
    for (const auto& matrix : matrices) {
        totalSize += matrix.N();
    }
    
    // Set up the merged matrix with random build mode (required for setrowsize)
    mergedMatrix.setSize(totalSize, totalSize);
    mergedMatrix.setBuildMode(WWMatrix::random);
    
    // Helper struct to find which original matrix a row belongs to
    struct MatrixLocation {
        int matrixIdx;
        int localRowIdx;
        int rowOffset;
    };
    
    auto findMatrixForRow = [&matrices](int globalRowIdx) -> MatrixLocation {
        int matrixIdx = 0;
        int localRowIdx = globalRowIdx;
        int rowOffset = 0;
        
        for (const auto& matrix : matrices) {
            if (localRowIdx < static_cast<int>(matrix.N())) {
                return {matrixIdx, localRowIdx, rowOffset};
            }
            localRowIdx -= matrix.N();
            rowOffset += matrix.N();
            matrixIdx++;
        }
        // Should never reach here if globalRowIdx is valid
        return {-1, -1, -1};
    };
    
    // Phase 1: Set row sizes
    for (int rowIdx = 0; rowIdx < totalSize; ++rowIdx) {
        auto loc = findMatrixForRow(rowIdx);
        if (loc.matrixIdx >= 0 && loc.matrixIdx < static_cast<int>(matrices.size())) {
            int nnz = countRowEntries(matrices[loc.matrixIdx][loc.localRowIdx]);
            mergedMatrix.setrowsize(rowIdx, nnz);
        }
    }
    mergedMatrix.endrowsizes();
    
    // Phase 2: Set column indices
    for (int rowIdx = 0; rowIdx < totalSize; ++rowIdx) {
        auto loc = findMatrixForRow(rowIdx);
        if (loc.matrixIdx >= 0 && loc.matrixIdx < static_cast<int>(matrices.size())) {
            const auto& origRow = matrices[loc.matrixIdx][loc.localRowIdx];
            for (auto colIt = origRow.begin(); colIt != origRow.end(); ++colIt) {
                mergedMatrix.addindex(rowIdx, loc.rowOffset + colIt.index());
            }
        }
    }
    mergedMatrix.endindices();
    
    // Phase 3: Copy values
    int rowOffset = 0;
    for (const auto& matrix : matrices) {
        for (int i = 0; i < static_cast<int>(matrix.N()); ++i) {
            for (auto colIt = matrix[i].begin(); colIt != matrix[i].end(); ++colIt) {
                mergedMatrix[rowOffset + i][rowOffset + colIt.index()] = *colIt;
            }
        }
        rowOffset += matrix.N();
    }
}

// Class to handle well matrix merging
class WellMatrixMerger
{
public:
    WellMatrixMerger() = default;
    
    // Add a well to the merger
    void addWell(const WRMatrix& B, const RWMatrix& C, const WWMatrix& D, 
                 const std::vector<int>& cellIndices, int wellIndex, const std::string& wellName)
    {
        B_matrices_.push_back(B);
        C_matrices_.push_back(C);
        D_matrices_.push_back(D);
        wellIndices_.push_back(wellIndex);
        wellNames_.push_back(wellName);
        
        // Map cell indices to well indices
        for (const auto& cellIdx : cellIndices) {
            cellToWellMap_[cellIdx] = wellIndex;
        }
    }
    
    // Finalize the merger by creating the merged matrices
    void finalize()
    {
        // Merge B matrices (well-to-reservoir)
        mergeWRMatrices(B_matrices_, mergedB_);
        
        // Merge C matrices (reservoir-to-well)
        mergeRWMatrices(C_matrices_, mergedC_);
        
        // Create diagonal block matrix for D matrices (well-to-well)
        createDiagonalBlockMatrix(D_matrices_, mergedD_);
    }
    
    // Get the merged matrices
    const WRMatrix& getMergedB() const { return mergedB_; }
    const RWMatrix& getMergedC() const { return mergedC_; }
    const WWMatrix& getMergedD() const { return mergedD_; }
    
    // Get well index for a cell
    int getWellIndexForCell(int cellIndex) const
    {
        auto it = cellToWellMap_.find(cellIndex);
        if (it != cellToWellMap_.end()) {
            return it->second;
        }
        return -1; // Not found
    }
    
    // Get well name by index
    const std::string& getWellName(int wellIndex) const
    {
        static const std::string notFound = "Unknown";
        for (size_t i = 0; i < wellIndices_.size(); ++i) {
            if (wellIndices_[i] == wellIndex) {
                return wellNames_[i];
            }
        }
        return notFound;
    }
    
    // Create system matrix
    SystemMatrix createSystemMatrix(const RRMatrix& A) const
    {
        // Create a system matrix with the correct structure
        SystemMatrix S;
        
        // Set the submatrices
        S[_0][_0] = A;
        S[_0][_1] = mergedC_;
        S[_1][_0] = mergedB_;
        S[_1][_1] = mergedD_;
        
        return S;
    }
    
private:
    std::vector<WRMatrix> B_matrices_;
    std::vector<RWMatrix> C_matrices_;
    std::vector<WWMatrix> D_matrices_;
    std::vector<int> wellIndices_;
    std::vector<std::string> wellNames_;
    std::map<int, int> cellToWellMap_;
    
    WRMatrix mergedB_;
    RWMatrix mergedC_;
    WWMatrix mergedD_;
};

// Helper for creating field matrices with a pattern
template<int ROWS, int COLS>
Dune::FieldMatrix<double, ROWS, COLS> makeBlock(const std::function<double(int, int)>& valueFn) {
    Dune::FieldMatrix<double, ROWS, COLS> block;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            block[i][j] = valueFn(i, j);
        }
    }
    return block;
}

// Helper for creating diagonal blocks with constant off-diagonal values
template<int SIZE>
Dune::FieldMatrix<double, SIZE, SIZE> makeDiagonalBlock(double diagValue, double offDiagValue = 0.0) {
    return makeBlock<SIZE, SIZE>([=](int i, int j) {
        return (i == j) ? diagValue : offDiagValue;
    });
}

// Helper for creating scaled blocks
template<int ROWS, int COLS>
Dune::FieldMatrix<double, ROWS, COLS> makeScaledBlock(double scale) {
    return makeBlock<ROWS, COLS>([=](int i, int j) {
        return scale * (i + 1) * (j + 1);
    });
}

// Create a tridiagonal matrix
void createTridiagonalMatrix(RRMatrix& matrix, int numRows, double diagValue, double offDiagValue) {
    matrix.setSize(numRows, numRows);
    matrix.setBuildMode(RRMatrix::random);
    
    // Phase 1: Set row sizes
    for (int i = 0; i < numRows; ++i) {
        int nnz = 1; // Diagonal element
        if (i > 0) nnz++; // Sub-diagonal
        if (i < numRows - 1) nnz++; // Super-diagonal
        matrix.setrowsize(i, nnz);
    }
    matrix.endrowsizes();
    
    // Phase 2: Set column indices
    for (int i = 0; i < numRows; ++i) {
        if (i > 0) {
            matrix.addindex(i, i-1); // Sub-diagonal
        }
        matrix.addindex(i, i); // Diagonal
        if (i < numRows - 1) {
            matrix.addindex(i, i+1); // Super-diagonal
        }
    }
    matrix.endindices();
    
    // Phase 3: Set values
    for (int i = 0; i < numRows; ++i) {
        if (i > 0) {
            matrix[i][i-1] = offDiagValue; // Sub-diagonal
        }
        matrix[i][i] = diagValue; // Diagonal
        if (i < numRows - 1) {
            matrix[i][i+1] = offDiagValue; // Super-diagonal
        }
    }
}

int main()
{
    std::cout << "=== Matrix Merger Example ===" << std::endl << std::endl;
    
    // Create matrices for a simple reservoir system
    RRMatrix A;
    
    // Create matrices for two wells
    WRMatrix B1, B2;
    RWMatrix C1, C2;
    WWMatrix D1, D2;
    
    // Cell indices for the wells
    std::vector<int> cells1 = {0, 1, 2}; // Cells for well 1
    std::vector<int> cells2 = {3, 4};    // Cells for well 2

    
    // Create reservoir matrix A (diagonal blocks)
    createTridiagonalMatrix(A, 5, 2.0, -1.0);

    // Create B1 matrix (well 1 to reservoir)
    auto B1_block = makeScaledBlock<numWellDofs, numResDofs>(0.1);
    createMatrix(B1, 1, 5, cells1, B1_block);
    
    // Create C1 matrix (reservoir to well 1)
    auto C1_block = makeScaledBlock<numResDofs, numWellDofs>(0.2);
    createMatrix(C1, 5, 1, {0}, C1_block);
    
    // Create D1 matrix (well 1 to well 1)
    auto D1_block = makeDiagonalBlock<numWellDofs>(1.5, 0.1);
    createMatrix(D1, 1, 1, {0}, D1_block);
    
    // Create B2 matrix (well 2 to reservoir)
    auto B2_block = makeScaledBlock<numWellDofs, numResDofs>(0.3);
    createMatrix(B2, 1, 5, cells2, B2_block);
    
    // Create C2 matrix (reservoir to well 2)
    auto C2_block = makeScaledBlock<numResDofs, numWellDofs>(0.4);
    createMatrix(C2, 5, 1, {0}, C2_block);
    
    // Create D2 matrix (well 2 to well 2)
    auto D2_block = makeDiagonalBlock<numWellDofs>(2.5, 0.2);
    createMatrix(D2, 1, 1, {0}, D2_block);
    
    // Print individual well matrices
    std::cout << "Individual well matrices:" << std::endl;
    printMatrix("A", A);
    printMatrix("B1", B1);
    printMatrix("C1", C1);
    printMatrix("D1", D1);
    printMatrix("B2", B2);
    printMatrix("C2", C2);
    printMatrix("D2", D2);
    
    // Create the well matrix merger
    WellMatrixMerger merger;
    
    // Add wells to the merger
    std::cout << "\nAdding wells to the merger..." << std::endl;
    merger.addWell(B1, C1, D1, cells1, 0, "Well1");
    merger.addWell(B2, C2, D2, cells2, 1, "Well2");
    
    // Finalize the merger
    std::cout << "Finalizing the merger..." << std::endl;
    merger.finalize();
    
    // Get the merged matrices
    const auto& mergedB = merger.getMergedB();
    const auto& mergedC = merger.getMergedC();
    const auto& mergedD = merger.getMergedD();
    
    // Print merged matrices
    std::cout << "\nMerged matrices:" << std::endl;
    printMatrix("mergedB", mergedB);
    printMatrix("mergedC", mergedC);
    printMatrix("mergedD", mergedD);
    
    // Check well indices for cells
    std::cout << "\nWell indices for cells:" << std::endl;
    for (const auto& cell : cells1) {
        std::cout << "Cell " << cell << " belongs to well " << merger.getWellIndexForCell(cell) 
                  << " (" << merger.getWellName(merger.getWellIndexForCell(cell)) << ")" << std::endl;
    }
    for (const auto& cell : cells2) {
        std::cout << "Cell " << cell << " belongs to well " << merger.getWellIndexForCell(cell)
                  << " (" << merger.getWellName(merger.getWellIndexForCell(cell)) << ")" << std::endl;
    }
    std::cout << "Cell 10 belongs to well " << merger.getWellIndexForCell(10) << " (not perforated)" << std::endl;
    
    // Create system matrix using the merged matrices
    std::cout << "\nCreating a system matrix using the merged matrices..." << std::endl;
    
    // Create a copy of the merged matrices to use in the system matrix
    RRMatrix A_copy = A;
    RWMatrix C_copy = mergedC;
    WRMatrix B_copy = mergedB;
    WWMatrix D_copy = mergedD;
    
    // Create the system matrix with the copied matrices
    SystemMatrix S;
    S[_0][_0] = A_copy;
    S[_0][_1] = C_copy;
    S[_1][_0] = B_copy;
    S[_1][_1] = D_copy;
    
    std::cout << "System matrix structure:" << std::endl;
    std::cout << "  A: " << A_copy.N() << " x " << A_copy.M() << std::endl;
    std::cout << "  C: " << C_copy.N() << " x " << C_copy.M() << std::endl;
    std::cout << "  B: " << B_copy.N() << " x " << B_copy.M() << std::endl;
    std::cout << "  D: " << D_copy.N() << " x " << D_copy.M() << std::endl;
    
    // SOLVER EXAMPLE
    std::cout << "\n=== Solver Example ===" << std::endl;
    
    // Create vectors for the system
    RVector x_r(5);
    WVector x_w(2);
    RVector r_res(5);
    WVector w_res(2);
    
    // Initialize vectors
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < numResDofs; ++j) {
            x_r[i][j] = 0.1 * (i + 1) * (j + 1);
            r_res[i][j] = 1.0;
        }
    }
    
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < numWellDofs; ++j) {
            x_w[i][j] = 0.2 * (i + 1) * (j + 1);
            w_res[i][j] = 1.0;
        }
    }
    
    // Create system vectors
    SystemVector x {x_r, x_w};
    SystemVector residual {r_res, w_res};
    
    printVector("Initial x_r", x[_0]);
    printVector("Initial x_w", x[_1]);
    printVector("Residual r_res", residual[_0]);
    printVector("Residual w_res", residual[_1]);
    
    // Create linear operator and preconditioner
    const Dune::MatrixAdapter<SystemMatrix, SystemVector, SystemVector> S_linop(S);
    TailoredPrecondDiag precond(S);
    
    // Set solver parameters
    double linsolve_tol = 1e-6;  // Less strict tolerance
    int max_iter = 20;           // Limit iterations
    int verbosity = 1;           // Reduce output verbosity
    
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
        solver.apply(x, residual, result);
        
        // Print results
        std::cout << "\nSolver results:" << std::endl;
        std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << std::endl;
        std::cout << "  Iterations: " << result.iterations << std::endl;
        std::cout << "  Reduction: " << result.reduction << std::endl;
        std::cout << "  Elapsed time: " << result.elapsed << " seconds" << std::endl;
        
        printVector("Solution x_r", x[_0]);
        printVector("Solution x_w", x[_1]);
    }
    catch (const std::exception& e) {
        std::cerr << "Error solving system: " << e.what() << std::endl;
    }
    
    std::cout << "\nMatrix merger example completed successfully!" << std::endl;
    
    return 0;
}


