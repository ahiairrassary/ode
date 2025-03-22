#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include <cstdint>
#include <functional>
#include <vector>

class AbstractOdeSolver {
public:
  virtual ~AbstractOdeSolver() = default;

  virtual float advance(size_t n) const = 0;

  void solve(float t0, float u0, float tFinal, uint32_t num) {
    m_t0 = t0;
    m_u0 = u0;
    m_tFinal = tFinal;
    m_num = num;
    m_dt = tFinal / num;

    m_t = std::vector<float>(num + 1, 0.0f);
    m_u = std::vector<float>(num + 1, 0.0f);

    m_t[0] = m_t0;
    m_u[0] = m_u0;

    for (size_t n = 0; n < num; ++n) {
      m_t[n + 1] = m_t[n] + m_dt;
      m_u[n + 1] = advance(n);
    }
  }

  std::vector<float> getT() const {
    return m_t;
  };

  std::vector<float> getU() const {
    return m_u;
  };

protected:
  std::function<float(float t, float u)> m_func;

  float m_t0;
  float m_u0;
  float m_tFinal;
  float m_num;
  float m_dt;

  std::vector<float> m_t;
  std::vector<float> m_u;
};

class ForwardEuler_v0 : public AbstractOdeSolver {
public:
  ForwardEuler_v0(const std::function<float(float t, float u)> &func) {
    m_func = func;
  }

  float advance(size_t n) const override {
    return m_u[n] + m_dt * m_func(m_t[n], m_u[n]);
  }
};

#endif // ODE_SOLVER_HPP
