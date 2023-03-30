#include "interface/LinearAlgebra/Vector.h"

namespace plasmatic {

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
Vector::Vector(Integer global_size) {
    PetscErrorCode ierr = VecCreate(PETSC_COMM_WORLD, &_data);

    if (ierr != 0) {
        std::abort();
    }

    ierr = VecSetSizes(_data, PETSC_DECIDE, global_size);
    if (ierr != 0) {
        std::abort();
    }

    ierr = VecSetFromOptions(_data);
    if (ierr != 0) {
        std::abort();
    }

    ierr = VecSet(_data, 0.0);
    if (ierr != 0) {
        std::abort();
    }
}

Integer Vector::Size() const {
    Integer size = 0;
    const PetscErrorCode ierr = VecGetSize(_data, &size);
    if (ierr != 0) {
        std::abort();
    }

    return size;
}

void Vector::AddValue(Integer pos, Float value) {
    const PetscErrorCode ierr = VecSetValues(_data, 1, &pos, &value, ADD_VALUES);
    if (ierr != 0) {
        std::abort();
    }
}

void Vector::SetValue(Integer pos, Float value) {
    const PetscErrorCode ierr = VecSetValues(_data, 1, &pos, &value, INSERT_VALUES);
    if (ierr != 0) {
        std::abort();
    }
}

void Vector::Assemble() {
    PetscErrorCode ierr = VecAssemblyBegin(_data);
    if (ierr != 0) {
        std::abort();
    }

    ierr = VecAssemblyEnd(_data);
    if (ierr != 0) {
        std::abort();
    }
}

Float Vector::GetValue(Integer pos) {
    Float value = std::numeric_limits<Float>::quiet_NaN();
    const PetscErrorCode ierr = VecGetValues(_data, 1, &pos, &value);
    if (ierr != 0) {
        std::abort();
    }

    return value;
}

} // namespace plasmatic
