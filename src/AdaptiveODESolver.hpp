#ifndef ADAPTIVE_ODE_SOLVER_HPP
#define ADAPTIVE_ODE_SOLVER_HPP

#include <functional>
#include <vector>

#include "ArrayLib.hpp"

template <size_t N, typename T = double>
class AbstractAdaptiveOdeSolver {
public:
  AbstractAdaptiveOdeSolver(const std::function<std::array<T, N>(T t, const std::array<T, N> &)> &func) {
    this->m_func = func;
  }

  virtual ~AbstractAdaptiveOdeSolver() = default;

  virtual T order() const = 0;
  virtual std::array<T, N> advance(T *pLocalError) const = 0;

  T firstStepSize(T t0, std::array<T, N> u0) const {
    const T DT = 0.1 / norm(m_func(t0, u0));

    return DT;
  }

  T newStepSize(T dt, T localError) const {
    if (std::isnan(localError) || std::isinf(localError)) {
      return m_minDt;

    } else {
      const T NEW_DT = (m_eta * std::pow(m_tol / localError, 1.0 / (order() + 1.0))) * dt;

      return std::clamp(NEW_DT, m_minDt, m_maxDt);
    }
  }

  void solve(T t0, T tFinal, std::array<T, N> u0, T tol = 1e-3, T minDt = 1e-5, T maxDt = 1e12) {
    m_t0 = t0;
    m_tFinal = tFinal;
    m_u0 = u0; // Initial condition

    m_tol = tol;
    m_minDt = minDt;
    m_maxDt = maxDt;

    m_t.reserve(1000); // TODO to optimize
    m_u.reserve(1000); // TODO to optimize

    m_t.resize(1);
    m_u.resize(1);

    m_t[0] = m_t0;
    m_u[0] = m_u0;

    T localT = m_t0;
    m_dt = firstStepSize(m_t0, m_u0);

    std::array<T, N> uNew;
    T localError;

    while (localT < m_tFinal) {
      uNew = advance(&localError);

      if (localError < m_tol || m_dt < m_minDt) {
        localT += m_dt;
        m_t.emplace_back(localT);
        m_u.emplace_back(uNew);

        m_dt = newStepSize(m_dt, localError);
        m_dt = std::min(m_dt, std::min(m_tFinal - localT, m_maxDt));

      } else {
        m_dt = newStepSize(m_dt, localError);
      }
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

  T m_eta = 0.9;

  T m_t0;
  T m_tFinal;

  T m_tol;
  T m_minDt;
  T m_maxDt;

  T m_dt;

  std::array<T, N> m_u0;

  std::vector<T> m_t;
  std::vector<std::array<T, N>> m_u;
};

template <size_t N, typename T = double>
class RKF45 : public AbstractAdaptiveOdeSolver<N, T> {
  // Runge–Kutta–Fehlberg 4(5)
public:
  using AbstractAdaptiveOdeSolver<N, T>::AbstractAdaptiveOdeSolver; // Constructor

  virtual T order() const override {
    return 4;
  }

  std::array<T, N> advance(T *pLocalError) const override {
    const T c2 = 1.0 / 4.0;
    const T a21 = 1.0 / 4.0;
    const T c3 = 3.0 / 8.0;
    const T a31 = 3.0 / 32.0;
    const T a32 = 9.0 / 32.0;
    const T c4 = 12.0 / 13.0;
    const T a41 = 1932.0 / 2197.0;
    const T a42 = -7200.0 / 2197.0;
    const T a43 = 7296.0 / 2197.0;
    const T c5 = 1.0;
    const T a51 = 439.0 / 216.0;
    const T a52 = -8.0;
    const T a53 = 3680.0 / 513.0;
    const T a54 = -845.0 / 4104.0;
    const T c6 = 1.0 / 2.0;
    const T a61 = -8.0 / 27.0;
    const T a62 = 2.0;
    const T a63 = -3544.0 / 2565.0;
    const T a64 = 1859.0 / 4104.0;
    const T a65 = -11.0 / 40.0;
    const T b1 = 25.0 / 216.0;
    // const T b2 = 0.0;
    const T b3 = 1408.0 / 2565.0;
    const T b4 = 2197.0 / 4104.0;
    const T b5 = -1.0 / 5.0;
    // const T b6 = 0.0;
    const T bh1 = 16.0 / 135.0;
    // const T bh2 = 0.0;
    const T bh3 = 6656.0 / 12825.0;
    const T bh4 = 28561.0 / 56430.0;
    const T bh5 = -9.0 / 50.0;
    const T bh6 = 2.0 / 55.0;

    const auto lastT = this->m_t.back();
    const auto lastU = this->m_u.back();

    const auto k1 = this->m_func(lastT, lastU);
    const auto k2 = this->m_func(lastT + c2 * this->m_dt, lastU + this->m_dt * (a21 * k1));
    const auto k3 = this->m_func(lastT + c3 * this->m_dt, lastU + this->m_dt * (a31 * k1 + a32 * k2));
    const auto k4 = this->m_func(lastT + c4 * this->m_dt, lastU + this->m_dt * (a41 * k1 + a42 * k2 + a43 * k3));
    const auto k5 =
        this->m_func(lastT + c5 * this->m_dt, lastU + this->m_dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4));
    const auto k6 = this->m_func(lastT + c6 * this->m_dt,
                                 lastU + this->m_dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5));

    const std::array<T, N> low = this->m_dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5);
    const std::array<T, N> high = this->m_dt * (bh1 * k1 + bh3 * k3 + bh4 * k4 + bh5 * k5 + bh6 * k6);

    std::array<T, N> unew = lastU + low;
    *pLocalError = norm(high - low);

    return unew;
  }
};

#endif // ADAPTIVE_ODE_SOLVER_HPP
