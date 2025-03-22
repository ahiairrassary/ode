#include <iostream>

#include "OdeSolver.hpp"

int main() {
  std::cout << "Hello ODE" << std::endl;

  auto logisticProblem = [](float t, float u) {
    (void)t;

    const float ALPHA = 0.2f;
    const float R = 1.0f;

    return ALPHA * u * (1 - u / R);
  };

  ForwardEuler_v0 FESolver(logisticProblem);
  FESolver.solve(0.0f, 0.1f, 40.0f, 400);

  const auto T = FESolver.getT();
  const auto U = FESolver.getU();

  for (size_t i = 0; i < T.size(); ++i) {
    std::cout << i << ": (" << T[i] << " ; " << U[i] << ")" << std::endl;
  }

  return EXIT_SUCCESS;
}
