#pragma once

#include "Utility/Utility.h"

#include <petscvec.h>

namespace plasmatic {

class Vector {
  public:
    Vector(Integer global_size);

    Integer Size() const;

    void AddValue(Integer pos, Float value);

    void SetValue(Integer pos, Float value);

    void Assemble();

    Float GetValue(Integer pos);

  private:
    Vec _data;
};

} // namespace plasmatic
