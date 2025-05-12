# Safe-Path-Planning with the Eikonal PDE

This project was built and tested on MATLAB R2024b.

MATLAB demo that converts **‖∇T‖ = c(x)** into collision-aware paths.

* **Maps** – three worlds  
  1. Single rectangle
  2. Four circles
  3. DFS-generated maze

* **Cost field** – `c(x) = 1 + α * exp(-β * dist(x, obs))`;  
  grid cells inside obstacles are assigned `1e5`.

* **Solver** – custom Fast-Marching implementation  
  (min-heap priority queue, 8-connected neighbors).

