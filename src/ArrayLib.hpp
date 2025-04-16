#ifndef ARRAY_LIB_HPP
#define ARRAY_LIB_HPP

#include <array>
#include <cmath>
#include <cstdint>

template <size_t N, typename T = double>
constexpr std::array<T, N> operator+(std::array<T, N> lhs, const std::array<T, N> &rhs) {
  for (size_t i = 0; i < N; ++i) {
    lhs[i] += rhs[i];
  }
  return lhs;
}

template <size_t N, typename T = double>
constexpr std::array<T, N> operator+(std::array<T, N> lhs, T rhs) {
  for (size_t i = 0; i < N; ++i) {
    lhs[i] += rhs;
  }
  return lhs;
}

template <size_t N, typename T = double>
constexpr std::array<T, N> operator-(std::array<T, N> lhs, const std::array<T, N> &rhs) {
  for (size_t i = 0; i < N; ++i) {
    lhs[i] -= rhs[i];
  }
  return lhs;
}

template <size_t N, typename T = double>
constexpr std::array<T, N> operator*(std::array<T, N> lhs, T rhs) {
  for (size_t i = 0; i < N; ++i) {
    lhs[i] *= rhs;
  }
  return lhs;
}

template <size_t N, typename T = double>
constexpr std::array<T, N> operator*(T lhs, std::array<T, N> rhs) {
  return rhs * lhs;
}

template <size_t N, typename T = double>
constexpr T norm(std::array<T, N> a) {
  T squareNorm = 0.0;

  for (size_t i = 0; i < N; ++i) {
    squareNorm += a[i] * a[i];
  }

  return std::sqrt(squareNorm);
}

#endif // ARRAY_LIB_HPP
