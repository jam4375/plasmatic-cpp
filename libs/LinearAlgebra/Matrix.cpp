#include "interface/LinearAlgebra/Matrix.h"

#include <petscksp.h>

namespace plasmatic {

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
Matrix::Matrix(Integer global_rows, Integer global_cols) {
    PetscErrorCode ierr = MatCreate(PETSC_COMM_WORLD, &_data);

    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    ierr = MatSetSizes(_data, PETSC_DECIDE, PETSC_DECIDE, global_rows, global_cols);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    ierr = MatSetFromOptions(_data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    ierr = MatSetUp(_data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
Matrix::Matrix(const Matrix &other) {
    const PetscErrorCode ierr = MatDuplicate(other._data, MAT_COPY_VALUES, &_data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

Matrix::~Matrix() {
    const PetscErrorCode ierr = MatDestroy(&_data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

Integer Matrix::Rows() const {
    Integer rows = 0;
    Integer cols = 0;
    const PetscErrorCode ierr = MatGetSize(_data, &rows, &cols);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    return rows;
}

Integer Matrix::Cols() const {
    Integer rows = 0;
    Integer cols = 0;
    const PetscErrorCode ierr = MatGetSize(_data, &rows, &cols);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    return cols;
}

void Matrix::AddValue(Integer row, Integer col, Float value) {
    const PetscErrorCode ierr = MatSetValues(_data, 1, &row, 1, &col, &value, ADD_VALUES);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

void Matrix::SetValue(Integer row, Integer col, Float value) {
    const PetscErrorCode ierr = MatSetValues(_data, 1, &row, 1, &col, &value, INSERT_VALUES);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

void Matrix::Assemble() {
    PetscErrorCode ierr = MatAssemblyBegin(_data, MAT_FINAL_ASSEMBLY);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    ierr = MatAssemblyEnd(_data, MAT_FINAL_ASSEMBLY);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

Float Matrix::GetValue(Integer row, Integer col) {
    Float value = std::numeric_limits<Float>::quiet_NaN();
    const PetscErrorCode ierr = MatGetValue(_data, row, col, &value);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    return value;
}

Matrix Matrix::operator+(const Matrix &other) {
    Matrix result(*this);

    const PetscErrorCode ierr = MatAXPY(result._data, 1.0, other._data, SAME_NONZERO_PATTERN);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
    return result;
}

Matrix Matrix::operator-(const Matrix &other) {
    Matrix result(*this);

    const PetscErrorCode ierr = MatAXPY(result._data, -1.0, other._data, SAME_NONZERO_PATTERN);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
    return result;
}

Matrix &Matrix::operator+=(const Matrix &other) {
    const PetscErrorCode ierr = MatAXPY(this->_data, 1.0, other._data, SAME_NONZERO_PATTERN);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
    return *this;
}

Matrix &Matrix::operator-=(const Matrix &other) {
    const PetscErrorCode ierr = MatAXPY(this->_data, -1.0, other._data, SAME_NONZERO_PATTERN);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
    return *this;
}

Vector Matrix::operator*(const Vector &other) {
    Vector result(other);

    const PetscErrorCode ierr = MatMult(this->_data, other._data, result._data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
    return result;
}

Vector Matrix::Solve(const Vector &other) {
    KSP ksp = nullptr;

    PetscErrorCode ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    ierr = KSPSetOperators(ksp, this->_data, this->_data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    PC preconditioner = nullptr;
    ierr = KSPGetPC(ksp, &preconditioner);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
    ierr = PCSetType(preconditioner, PCJACOBI);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
    constexpr auto tol = 1.0e-10;
    ierr = KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    ierr = KSPSetFromOptions(ksp);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    Vector result(other);
    ierr = KSPSolve(ksp, other._data, result._data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    ierr = KSPDestroy(&ksp);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    return result;
}

void Matrix::SetDirichletBC(Integer row_col, const Vector &x, const Vector &b) {
    const PetscErrorCode ierr = MatZeroRowsColumns(_data, 1, &row_col, 1.0, x._data, b._data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

} // namespace plasmatic
