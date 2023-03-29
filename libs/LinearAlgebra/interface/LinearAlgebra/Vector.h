#pragma once

#include "Utility/Utility.h"

#include <petscvec.h>

namespace plasmatic {

class Vector {
  public:
    Vector(Integer global_size);

    Integer Size() const;

  private:
    Vec _data;
};

} // namespace plasmatic
