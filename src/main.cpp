#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "AdaptiveODESolver.hpp"
#include "OdeSolver.hpp"

#define ASSERT_DOUBLE(a, b) (assert(std::abs(a - b) < 1e-15))

int main() {
  std::cout << "Hello ODE" << std::endl;
  std::cout << std::setprecision(16);

  auto logisticProblem = [](double t, const std::array<double, 1> &u) {
    (void)t;

    // Inputs
    const double ALPHA = 0.2;
    const double R = 1.0;

    std::array<double, 1> tmp;

    for (size_t i = 0; i < tmp.size(); ++i) {
      tmp[i] = ALPHA * u[i] * (1.0 - u[i] / R);
    }

    return tmp;
  };

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

  {
    ForwardEuler<1, double> FESolver(logisticProblem);
    FESolver.solve(0.0, 40.0, 400, {0.1});

    /*const auto T = FESolver.getT();
    const auto U = FESolver.getU();

    for (size_t i = 0; i < T.size(); ++i) {
      std::cout << i << ": (" << T[i] << " ; " << U[i][0] << ")" << std::endl;
    }*/
  }

  {
    RungeKutta4<2, double> FESolver(pendulumProblem);
    FESolver.solve(0.0, 10.0, 1000.0, {M_PI / 4.0, 0.0});

    /*const auto T = FESolver.getT();
    const auto U = FESolver.getU();

    for (size_t i = 0; i < T.size(); ++i) {
      std::cout << i << ": t=" << T[i] << ": (" << U[i][0] << " : " << U[i][1] << ")" << std::endl;
    }*/
  }

  {
    RKF45<2, double> RKF45Solver(pendulumProblem);
    RKF45Solver.solve(0.0, 10.0, {M_PI / 4.0, 0.0}, 1e-6);

    const auto T = RKF45Solver.getT();
    const auto U = RKF45Solver.getU();

    assert(T.size() == 147);

    for (size_t i = 0; i < T.size(); ++i) {
      std::cout << std::setw(4) << i << ": t=" << std::setw(20) << T[i] << ": (" << std::setw(20) << U[i][0] << " : "
                << std::setw(20) << U[i][1] << ")" << std::endl;
    }

    ASSERT_DOUBLE(T[0], 0.0);
    ASSERT_DOUBLE(U[0][0], 0.7853981633974483);
    ASSERT_DOUBLE(U[0][1], 0.0);

    ASSERT_DOUBLE(T[1], 0.014416040391163046);
    ASSERT_DOUBLE(U[1][0], 0.7846774479842169);
    ASSERT_DOUBLE(U[1][1], -0.09997597008360871);

    ASSERT_DOUBLE(T[2], 0.08108974990631218);
    ASSERT_DOUBLE(U[2][0], 0.7626788001049662);
    ASSERT_DOUBLE(U[2][1], -0.5582019687198568);

    ASSERT_DOUBLE(T[145], 9.982343918355255);
    ASSERT_DOUBLE(U[145][0], 0.1726281618263992);
    ASSERT_DOUBLE(U[145][1], 2.335658033834428);

    ASSERT_DOUBLE(T[146], 10.0);
    ASSERT_DOUBLE(U[146][0], 0.2135834808194253);
    ASSERT_DOUBLE(U[146][1], 2.302412409767315);
  }

  return EXIT_SUCCESS;
}
