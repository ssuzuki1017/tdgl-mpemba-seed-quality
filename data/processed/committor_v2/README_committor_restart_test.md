# Committor-style restart test

This directory contains a finite-time committor-style restart analysis.

Parameters:

```text
N = 64
dt = 0.02
tmax = 300.0
preeq_steps = 500
a_f = 0.02
D_f = 0.009
D0 = 0.02
cluster_threshold = 20
min_consecutive = 3
n_configs_per_label = 12
n_restarts = 24
labels = (1.05, 1.1, 1.2, 1.5, 2.0, 3.0)
```

Definition:

```text
q_nuc(phi0) = fraction of independent post-quench restarts that nucleate by tmax
```

Outputs:

- `committor_initial_configs.csv`
- `committor_restarts.csv`
- `committor_label_summary.csv`
- `committor_geometry_correlations.csv`
- `figS7_committor_seed_geometry.png`
- `figS7_committor_seed_geometry.pdf`

Interpretation:

This is not an exact transition-path committor because the observation time is finite and the nucleation event is defined by an operational cluster threshold. It is nevertheless a direct restart test of whether initial seed geometry contains information about finite-time nucleation outcomes.
