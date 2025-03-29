#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include <cstdint>
#include <functional>
#include <vector>

template <size_t N, typename T = double>
class AbstractOdeSolver {
public:
  virtual ~AbstractOdeSolver() = default;

  virtual std::array<T, N> advance(size_t n) const = 0;

  void solve(T t0, T tFinal, uint32_t num, std::array<T, N> u0) {
    m_t0 = t0;
    m_tFinal = tFinal;
    m_num = num;
    m_dt = (tFinal - t0) / num;

    m_u0 = u0;

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
class ForwardEuler_v0 : public AbstractOdeSolver<N, T> {
public:
  ForwardEuler_v0(const std::function<std::array<T, N>(T t, const std::array<T, N> &)> &func) {
    this->m_func = func;
  }

  std::array<T, N> advance(size_t n) const override {
    const auto F_N = this->m_func(this->m_t[n], this->m_u[n]);

    std::array<T, N> tmp;

    for (size_t i = 0; i < tmp.size(); ++i) {
      tmp[i] = this->m_u[n][i] + this->m_dt * F_N[i];
    }

    return tmp;
  }
};

#endif // ODE_SOLVER_HPP
