#pragma once

#include "Utility/Utility.h"

#include "Vector.h"

#include <petscmat.h>

namespace plasmatic {

class Matrix {
  public:
    Matrix(Integer global_rows, Integer global_cols);

    Matrix(const Matrix &other);

    Integer Rows() const;

    Integer Cols() const;

    void AddValue(Integer row, Integer col, Float value);

    void SetValue(Integer row, Integer col, Float value);

    void Assemble();

    Float GetValue(Integer row, Integer col);

    Matrix operator+(const Matrix &other);

    Matrix operator-(const Matrix &other);

    Matrix &operator+=(const Matrix &other);

    Matrix &operator-=(const Matrix &other);

    Vector operator*(const Vector &other);

    Vector Solve(const Vector &other);

  private:
    Mat _data;
};

} // namespace plasmatic
