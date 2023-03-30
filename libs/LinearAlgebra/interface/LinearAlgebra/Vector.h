#pragma once

#include "Utility/Utility.h"

#include <petscvec.h>

namespace plasmatic {

class Vector {
  public:
    Vector(Integer global_size);

    Vector(Vector &other);

    Integer Size() const;

    void AddValue(Integer pos, Float value);

    void SetValue(Integer pos, Float value);

    void Assemble();

    Float GetValue(Integer pos);

    Vector operator+(const Vector &other);

    Vector operator-(const Vector &other);

    Vector &operator+=(const Vector &other);

    Vector &operator-=(const Vector &other);

  private:
    Vec _data;
};

} // namespace plasmatic
