#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include <functional>
#include <vector>

#include "ArrayLib.hpp"

// ****************************************

template <size_t N, typename T>
class AbstractOdeSolver {
public:
  AbstractOdeSolver(const std::function<std::array<T, N>(T t, const std::array<T, N> &)> &func) {
    this->m_func = func;
  }

  virtual ~AbstractOdeSolver() = default;

  // Compute u[n+1]
  virtual std::array<T, N> advance(size_t n) const = 0;

  void solve(T t0, T tFinal, uint32_t num, const std::array<T, N> &u0) {
    m_t0 = t0;
    m_tFinal = tFinal;
    m_num = num;

    m_u0 = u0; // Initial condition

    m_dt = (tFinal - t0) / num;

    m_t.resize(num + 1);
    m_u.resize(num + 1);

    m_t[0] = m_t0;
    m_u[0] = m_u0;

    for (size_t n = 0; n < num; ++n) {
      m_t[n + 1] = m_t[n] + m_dt;
      m_u[n + 1] = advance(n);
    }
  }

  std::vector<T> getT() const {
    return m_t;
  };

  std::vector<std::array<T, N>> getU() const {
    return m_u;
  };

protected:
  std::function<std::array<T, N>(T t, const std::array<T, N> &)> m_func;

  T m_t0;
  T m_tFinal;
  T m_num;

  T m_dt;

  std::array<T, N> m_u0;

  std::vector<T> m_t;
  std::vector<std::array<T, N>> m_u;
};

template <size_t N, typename T = double>
class ForwardEuler : public AbstractOdeSolver<N, T> {
public:
  using AbstractOdeSolver<N, T>::AbstractOdeSolver; // Constructor

  std::array<T, N> advance(size_t n) const override {
    return (this->m_u[n] + this->m_dt * this->m_func(this->m_t[n], this->m_u[n]));
  }
};

template <size_t N, typename T = double>
class RungeKutta4 : public AbstractOdeSolver<N, T> {
public:
  using AbstractOdeSolver<N, T>::AbstractOdeSolver; // Constructor

  std::array<T, N> advance(size_t n) const override {
    const T dt2 = this->m_dt / 2.0;

    const auto k1 = this->m_func(this->m_t[n], this->m_u[n]);
    const auto k2 = this->m_func(this->m_t[n] + dt2, this->m_u[n] + dt2 * k1);
    const auto k3 = this->m_func(this->m_t[n] + dt2, this->m_u[n] + dt2 * k2);
    const auto k4 = this->m_func(this->m_t[n] + this->m_dt, this->m_u[n] + this->m_dt * k3);

    return (this->m_u[n] + (this->m_dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4));
  }
};

#endif // ODE_SOLVER_HPP
