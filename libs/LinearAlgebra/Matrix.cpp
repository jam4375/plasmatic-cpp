#include "interface/LinearAlgebra/Matrix.h"

namespace plasmatic {

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
Matrix::Matrix(Integer global_rows, Integer global_cols) {
    PetscErrorCode ierr = MatCreate(PETSC_COMM_WORLD, &_data);

    if (ierr != 0) {
        std::abort();
    }

    ierr = MatSetSizes(_data, PETSC_DECIDE, PETSC_DECIDE, global_rows, global_cols);
    if (ierr != 0) {
        std::abort();
    }

    ierr = MatSetFromOptions(_data);
    if (ierr != 0) {
        std::abort();
    }

    ierr = MatSetUp(_data);
    if (ierr != 0) {
        std::abort();
    }
}

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
Matrix::Matrix(Matrix &other) {
    const PetscErrorCode ierr = MatDuplicate(other._data, MAT_COPY_VALUES, &_data);
    if (ierr != 0) {
        std::abort();
    }
}

Integer Matrix::Rows() const {
    Integer rows = 0;
    Integer cols = 0;
    const PetscErrorCode ierr = MatGetSize(_data, &rows, &cols);
    if (ierr != 0) {
        std::abort();
    }

    return rows;
}

Integer Matrix::Cols() const {
    Integer rows = 0;
    Integer cols = 0;
    const PetscErrorCode ierr = MatGetSize(_data, &rows, &cols);
    if (ierr != 0) {
        std::abort();
    }

    return cols;
}

void Matrix::AddValue(Integer row, Integer col, Float value) {
    const PetscErrorCode ierr = MatSetValues(_data, 1, &row, 1, &col, &value, ADD_VALUES);
    if (ierr != 0) {
        std::abort();
    }
}

void Matrix::SetValue(Integer row, Integer col, Float value) {
    const PetscErrorCode ierr = MatSetValues(_data, 1, &row, 1, &col, &value, INSERT_VALUES);
    if (ierr != 0) {
        std::abort();
    }
}

void Matrix::Assemble() {
    PetscErrorCode ierr = MatAssemblyBegin(_data, MAT_FINAL_ASSEMBLY);
    if (ierr != 0) {
        std::abort();
    }

    ierr = MatAssemblyEnd(_data, MAT_FINAL_ASSEMBLY);
    if (ierr != 0) {
        std::abort();
    }
}

Float Matrix::GetValue(Integer row, Integer col) {
    Float value = std::numeric_limits<Float>::quiet_NaN();
    const PetscErrorCode ierr = MatGetValue(_data, row, col, &value);
    if (ierr != 0) {
        std::abort();
    }

    return value;
}

Matrix Matrix::operator+(const Matrix &other) {
    Matrix result(*this);

    const PetscErrorCode ierr = MatAXPY(result._data, 1.0, other._data, SAME_NONZERO_PATTERN);
    if (ierr != 0) {
        std::abort();
    }
    return result;
}

Matrix Matrix::operator-(const Matrix &other) {
    Matrix result(*this);

    const PetscErrorCode ierr = MatAXPY(result._data, -1.0, other._data, SAME_NONZERO_PATTERN);
    if (ierr != 0) {
        std::abort();
    }
    return result;
}

Matrix &Matrix::operator+=(const Matrix &other) {
    const PetscErrorCode ierr = MatAXPY(this->_data, 1.0, other._data, SAME_NONZERO_PATTERN);
    if (ierr != 0) {
        std::abort();
    }
    return *this;
}

Matrix &Matrix::operator-=(const Matrix &other) {
    const PetscErrorCode ierr = MatAXPY(this->_data, -1.0, other._data, SAME_NONZERO_PATTERN);
    if (ierr != 0) {
        std::abort();
    }
    return *this;
}

} // namespace plasmatic
