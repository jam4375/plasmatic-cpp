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
}

Integer Vector::Size() const {
    Integer size = 0;
    const PetscErrorCode ierr = VecGetSize(_data, &size);
    if (ierr != 0) {
        std::abort();
    }

    return size;
}

} // namespace plasmatic
