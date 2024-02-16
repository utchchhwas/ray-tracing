#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <array>
#include <stdexcept>

// Matrix class
template <typename T, size_t Rows, size_t Cols>
class Matrix {

private:
    std::array<std::array<T, Cols>, Rows> data;

public:
    // Default constructor
    Matrix() : data{} {}

    // Copy constructor
    Matrix(const Matrix& other) {
        data = other.data;
    }

    // Constructor from initializer list
    Matrix(const std::initializer_list<std::initializer_list<T>>& values) {
        if (values.size() != Rows) {
            throw std::invalid_argument("Incorrect number of rows in initializer list.");
        }

        size_t i = 0;
        for (const auto& row : values) {
            if (row.size() != Cols) {
                throw std::invalid_argument("Incorrect number of columns in initializer list.");
            }

            size_t j = 0;
            for (const auto& value : row) {
                data[i][j] = value;
                ++j;
            }
            ++i;
        }
    }

    // Constructor from 1D initializer list
    Matrix(const std::initializer_list<T>& values) {
        if (values.size() != Rows * Cols) {
            throw std::invalid_argument("Incorrect number of elements in initializer list.");
        }

        size_t i = 0;
        size_t j = 0;
        for (const auto& value : values) {
            data[i][j] = value;
            ++j;
            if (j == Cols) {
                j = 0;
                ++i;
            }
        }
    }

    // Assignment operator
    Matrix& operator=(const Matrix& other) {
        if (this == &other) {
            return *this;
        }
        data = other.data;
        return *this;
    }

    // Get number of rows
    size_t getRows() const {
        return Rows;
    }

    // Get number of columns
    size_t getCols() const {
        return Cols;
    }

    // Overloaded [] operator for row access
    std::array<T, Cols>& operator[](size_t row) {
        if (row >= Rows) {
            throw std::out_of_range("Matrix row index out of range.");
        }
        return data[row];
    }

    // Overloaded [] operator for const row access
    const std::array<T, Cols>& operator[](size_t row) const {
        if (row >= Rows) {
            throw std::out_of_range("Matrix row index out of range.");
        }
        return data[row];
    }

    // Overloaded addition operator
    Matrix operator+(const Matrix& other) const {
        Matrix result;
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    // Overloaded multiplication operator
    template <size_t OtherCols>
    Matrix<T, Rows, OtherCols> operator*(const Matrix<T, Cols, OtherCols>& other) const {
        Matrix<T, Rows, OtherCols> result;
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < OtherCols; ++j) {
                T sum = 0;
                for (size_t k = 0; k < Cols; ++k) {
                    sum += (*this)[i][k] * other[k][j];
                }
                result[i][j] = sum;
            }
        }
        return result;
    }

    // Static method to get the identity matrix
    static Matrix<T, Rows, Cols> identity() {
        Matrix<T, Rows, Cols> result;
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                result.data[i][j] = (i == j) ? T(1) : T(0);
            }
        }
        return result;
    }

    // Method to calculate the determinant of a square matrix
    T determinant() const {
        static_assert(Rows == Cols, "Determinant is only defined for square matrices.");

        if (Rows == 1) {
            return data[0][0];
        } else if (Rows == 2) {
            return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        } else {
            T det = 0;
            for (size_t j = 0; j < Cols; ++j) {
                Matrix<T, Rows - 1, Cols - 1> submatrix;
                for (size_t i = 1; i < Rows; ++i) {
                    size_t subRow = 0;
                    for (size_t subCol = 0; subCol < Cols; ++subCol) {
                        if (subCol != j) {
                            submatrix[subRow][subCol] = data[i][subCol];
                            ++subRow;
                        }
                    }
                }
                det += (j % 2 == 0 ? 1 : -1) * data[0][j] * submatrix.determinant();
            }
            return det;
        }
    }
};

// Overloaded stream insertion operator for convenient output
template <typename T, size_t Rows, size_t Cols>
std::ostream& operator<<(std::ostream& os, const Matrix<T, Rows, Cols>& matrix) {
    for (size_t i = 0; i < Rows; ++i) {
        for (size_t j = 0; j < Cols; ++j) {
            os << matrix[i][j] << ' ';
        }
        if (i != Rows - 1) {
            os << '\n';
        }
    }
    return os;
}

#endif // MATRIX_HPP
