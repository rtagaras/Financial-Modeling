# Options pricing with Monte Carlo

This project demonstrates methods for pricing a variety of options. 
options_pricing.cpp calculates values for the following option types:

- European options (binomial lattice, geometric Brownian motion, jump-diffusion process)
- American options (binomial lattice)
- Fixed strike Asian options (geometric Brownian motion, jump-diffusion process)
- Barrier Options (geometric Brownian motion)
- Basket Options (geometric Brownian motion)
- Exchange Options (geometric Brownian motion)
- Bermudan Options (binomial lattice)

processing.py calculates confidence intervals for the option values computed by options_pricing.cpp

american_boundary.cpp prices an American option using simulated annealing to find an optimal exercise boundary. 
