#ifndef VEC_HPP
#define VEC_HPP

#include <iostream>
#include <array>
#include <stdexcept>
#include <cmath>

template <typename T, size_t Size>
class Vec {
private:
    std::array<T, Size> data;

public:
    Vec() : data{} {}

    Vec(const std::initializer_list<T>& values) {
        if (values.size() != Size) {
            throw std::invalid_argument("Incorrect number of elements in initializer list.");
        }

        size_t i = 0;
        for (const auto& value : values) {
            data[i] = value;
            ++i;
        }
    }

    T& operator[](size_t index) {
        if (index >= Size) {
            throw std::out_of_range("Vec index out of range.");
        }
        return data[index];
    }

    const T& operator[](size_t index) const {
        if (index >= Size) {
            throw std::out_of_range("Vec index out of range.");
        }
        return data[index];
    }

    // Vector negation
    Vec operator-() const {
        Vec<T, Size> result;
        for (size_t i = 0; i < Size; ++i) {
            result[i] = -data[i];
        }
        return result;
    }

    // Vector addition
    Vec operator+(const Vec<T, Size>& other) const {
        Vec<T, Size> result;
        for (size_t i = 0; i < Size; ++i) {
            result[i] = data[i] + other[i];
        }
        return result;
    }

    // Vector subtraction
    Vec operator-(const Vec<T, Size>& other) const {
        Vec<T, Size> result;
        for (size_t i = 0; i < Size; ++i) {
            result[i] = data[i] - other[i];
        }
        return result;
    }

    // Scalar multiplication (vector * scalar)
    Vec operator*(T scalar) const {
        Vec<T, Size> result;
        for (size_t i = 0; i < Size; ++i) {
            result[i] = data[i] * scalar;
        }
        return result;
    }

    // Scalar multiplication (scalar * vector)
    friend Vec operator*(T scalar, const Vec<T, Size>& vector) {
        Vec<T, Size> result;
        for (size_t i = 0; i < Size; ++i) {
            result[i] = vector[i] * scalar;
        }
        return result;
    }

    // Compute the magnitude (length) of the vector
    T magnitude() const {
        T sumSquares = 0;
        for (size_t i = 0; i < Size; ++i) {
            sumSquares += data[i] * data[i];
        }
        return std::sqrt(sumSquares);
    }

    // Normalize the vector
    void normalize() {
        T magnitude = this->magnitude();
        if (magnitude != 0) {
            T reciprocalMagnitude = 1 / magnitude;
            for (size_t i = 0; i < Size; ++i) {
                data[i] *= reciprocalMagnitude;
            }
        }
    }

    // Dot product
    T dot(const Vec<T, Size>& other) const {
        T result = 0;
        for (size_t i = 0; i < Size; ++i) {
            result += data[i] * other[i];
        }
        return result;
    }

    // Cross product (for 3D vectors only)
    Vec<T, 3> cross(const Vec<T, 3>& other) const {
        static_assert(Size == 3, "Cross product is only defined for 3D vectors.");

        return Vec<T, 3>{
            data[1] * other[2] - data[2] * other[1],
            data[2] * other[0] - data[0] * other[2],
            data[0] * other[1] - data[1] * other[0]
        };
    }

    // Rotate the vector using Rodrigues' rotation formula
    Vec<T, 3> rotate(const Vec<T, 3>& axis, double angleDegrees) const {
        // Convert angle to radians
        double angleRadians = angleDegrees * (M_PI / 180.0);

        // Calculate sin and cos of the angle
        double c = cos(angleRadians);
        double s = sin(angleRadians);
        double t = 1.0 - c;

        // Calculate the rotated vector
        Vec<T, 3> rotatedVector = *this * c + axis.cross(*this) * s + axis * axis.dot(*this) * t;

        return rotatedVector;
    }
};

template <typename T, size_t Size>
std::ostream& operator<<(std::ostream& os, const Vec<T, Size>& vector) {
    for (size_t i = 0; i < Size; ++i) {
        os << vector[i];
        if (i != Size - 1) {
            os << ' ';
        }
    }
    return os;
}

#endif // VEC_HPP
