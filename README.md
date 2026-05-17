# Travelling Salesman Problem — GA & Simulated Annealing

A Python implementation of two metaheuristic solutions to the **Travelling Salesman Problem (TSP)**: a **Genetic Algorithm (GA)** and **Simulated Annealing (SA)**. After each run, the best route found at every recorded generation is rendered step-by-step in Pygame, with a Tkinter slider for playback. Convergence data is exported as a CSV log and can be plotted with `plot.py`.

## Methods

### Genetic Algorithm
Evolves a population of route permutations across generations. Supports two parent selection strategies — **Rank Selection (RS)** and **Stochastic Universal Sampling (SUS)** — switchable via `PAR_SEL_METHOD`. Crossover uses the Edge Recombination Operator (60% probability per pair), and offspring are subject to one of 7 mutation operators at a 1–10% rate. The run halts early if improvement stagnates below a configurable threshold.

### Simulated Annealing
Optimises a single permutation by mutating it each iteration and accepting worse solutions with probability `P = e^(−ΔE / T)`. Temperature decays each generation, gradually reducing acceptance of worse solutions. The same 7 mutation operators are used as in the GA.

Switch between methods by setting `SOL_METHOD = "GA"` or `"SA"` in `TSP.py`.

## Project Structure

```
TSP.py                  # Main solver
plot.py                 # Convergence plotter
GA_succession_SUS.txt   # GA log — Stochastic Universal Sampling
GA_succession_RS.txt    # GA log — Rank Selection
SA_succession.txt       # SA convergence log
documentation.docx      # Full written documentation
```

## Requirements

```bash
pip install pygame matplotlib pandas numpy Pillow
```

## Usage

```bash
python TSP.py   # run the solver and visualize results
python plot.py  # plot convergence from a saved log file
```
