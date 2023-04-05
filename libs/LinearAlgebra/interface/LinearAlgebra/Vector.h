#pragma once

#include "Utility/Utility.h"

#include <petscvec.h>

namespace plasmatic {

class Vector {
  public:
    Vector(Integer global_size);

    Vector(const Vector &other);

    ~Vector();

    Integer Size() const;

    void AddValue(Integer pos, Float value);

    void SetValue(Integer pos, Float value);

    void Assemble();

    Float GetValue(Integer pos);

    Vector operator+(const Vector &other);

    Vector operator-(const Vector &other);

    Vector &operator+=(const Vector &other);

    Vector &operator-=(const Vector &other);

    friend class Matrix;

  private:
    Vec _data;
};

} // namespace plasmatic
