# Experiment test simulation 6.7.12

This experiment tests the performance of intrinsic MCMC (iMCMC).

Notes:
- N=(2000, 2000)
- phi = (0, 1, 10)
- `dataseed` = (1, 2, 3)
- `mcmcseed` = (1, 2, 3)
- omega ~ Normal(-3, 1.3) => p ~ Beta(1, 9) (weakly encourage sparsity)
- iMCMC will be fast, but we also hope it helps with mixing, and recovers the
  simulation truth well.

