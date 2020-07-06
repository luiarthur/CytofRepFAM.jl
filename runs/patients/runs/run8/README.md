# RUN 1

- Configs
    - Burn in: 10000
    - MCMC samples: 5000
    - Thinning: None
    - K = 20
    - L = 5
    - `W*` ~ Gamma(10, rate=1)
    - `p_i` ~ Beta(1, 99) => `omega_q` ~ Normal(-5, 1)
    - phi = (0, 1, 10, 100) (command line arg)
    - Temperatures = (1, 1.0003, 1.0006, 1.001)
    - Save all states
    - Patients Data: (001_d31, 007_d35, 010_d35)
    - removed markers:
        - 2, 4, 6, 9, 10, 20, 21, 23, 31
        - These markers are either mostly non-expressed or mostly expressed
    - Run with samples 1 and 2 only.

- Visualizations
    - [ ] TSNE
        - Color by cluster

- Comparisons
    - Not in paper
        - [ ] Flowsom
        - [ ] Mclust
