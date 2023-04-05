#include "interface/LinearAlgebra/Vector.h"

namespace plasmatic {

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
Vector::Vector(Integer global_size) {
    PetscErrorCode ierr = VecCreate(PETSC_COMM_WORLD, &_data);

    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    ierr = VecSetSizes(_data, PETSC_DECIDE, global_size);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    ierr = VecSetFromOptions(_data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    ierr = VecSet(_data, 0.0);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
Vector::Vector(const Vector &other) {
    PetscErrorCode ierr = VecDuplicate(other._data, &_data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    ierr = VecCopy(other._data, _data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

Vector::~Vector() {
    const PetscErrorCode ierr = VecDestroy(&_data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

Integer Vector::Size() const {
    Integer size = 0;
    const PetscErrorCode ierr = VecGetSize(_data, &size);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    return size;
}

void Vector::AddValue(Integer pos, Float value) {
    const PetscErrorCode ierr = VecSetValues(_data, 1, &pos, &value, ADD_VALUES);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

void Vector::SetValue(Integer pos, Float value) {
    const PetscErrorCode ierr = VecSetValues(_data, 1, &pos, &value, INSERT_VALUES);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

void Vector::Assemble() {
    PetscErrorCode ierr = VecAssemblyBegin(_data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    ierr = VecAssemblyEnd(_data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
}

Float Vector::GetValue(Integer pos) {
    Float value = std::numeric_limits<Float>::quiet_NaN();
    const PetscErrorCode ierr = VecGetValues(_data, 1, &pos, &value);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    return value;
}

Vector Vector::operator+(const Vector &other) {
    Vector result(this->Size());

    const PetscErrorCode ierr = VecAXPBYPCZ(result._data, 1.0, 1.0, 0.0, _data, other._data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
    return result;
}

Vector Vector::operator-(const Vector &other) {
    Vector result(this->Size());

    const PetscErrorCode ierr = VecAXPBYPCZ(result._data, 1.0, -1.0, 0.0, _data, other._data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
    return result;
}

Vector &Vector::operator+=(const Vector &other) {
    const PetscErrorCode ierr = VecAXPY(this->_data, 1.0, other._data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
    return *this;
}

Vector &Vector::operator-=(const Vector &other) {
    const PetscErrorCode ierr = VecAXPY(this->_data, -1.0, other._data);
    Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);
    return *this;
}

} // namespace plasmatic
