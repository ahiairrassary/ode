#include <cmath>
#include <iomanip>
#include <iostream>

#include "OdeSolver.hpp"

int main() {
  std::cout << "Hello ODE" << std::endl;
  std::cout << std::setprecision(15);

  {
    auto logisticProblem = [](double t, const std::array<double, 1> &u) {
      (void)t;

      const double ALPHA = 0.2;
      const double R = 1.0;

      std::array<double, 1> tmp;

      for (size_t i = 0; i < tmp.size(); ++i) {
        tmp[i] = ALPHA * u[i] * (1.0 - u[i] / R);
      }

      return tmp;
    };

    ForwardEuler<1, double> FESolver(logisticProblem);
    FESolver.solve(0.0, 40.0, 400.0, {0.1});

    /*const auto T = FESolver.getT();
    const auto U = FESolver.getU();

    for (size_t i = 0; i < T.size(); ++i) {
      std::cout << i << ": (" << T[i] << " ; " << U[i][0] << ")" << std::endl;
    }*/
  }

  auto pendulumProblem = [](double t, const std::array<double, 2> &u) {
    (void)t;

    const double L = 1.0;
    const double G = 9.81;

    const double theta = u[0];
    const double omega = u[1];

    std::array<double, 2> tmp;
    tmp[0] = omega;                    // dtheta
    tmp[1] = -G / L * std::sin(theta); // domega

    return tmp;
  };

  RungeKutta4<2, double> FESolver(pendulumProblem);
  FESolver.solve(0.0, 10.0, 1000.0, {M_PI / 4.0, 0.0});

  const auto T = FESolver.getT();
  const auto U = FESolver.getU();

  for (size_t i = 0; i < T.size(); ++i) {
    std::cout << i << ": t=" << T[i] << ": (" << U[i][0] << " : " << U[i][1] << ")" << std::endl;
  }

  return EXIT_SUCCESS;
}
