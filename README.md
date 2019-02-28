# Nonlinear pendulum solution with Crank-Nicolson finite difference method
The nonlinear pendulum problem θ"+sinθ=0 is solved by converting the equation into two equations dθ'/dt=-sinθ and dθ/dt=θ'. Then it is solved with finite difference using corrector-predictor and Crank-Nicolson method. The latter method fixes the nonlinearity of sinθ by using Taylor series up to the first term.
