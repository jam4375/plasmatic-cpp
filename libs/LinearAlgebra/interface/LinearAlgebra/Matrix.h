#pragma once

#include "Utility/Utility.h"

#include <petscmat.h>

namespace plasmatic {

class Matrix {
  public:
    Matrix(Integer global_rows, Integer global_cols);

    Matrix(Matrix &other);

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

  private:
    Mat _data;
};

} // namespace plasmatic
