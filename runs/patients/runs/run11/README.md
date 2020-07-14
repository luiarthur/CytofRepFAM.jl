# RUN 1

- Configs (same as run10, but with sb-ibp)
    - Burn in: 10000
    - MCMC samples: 5000
    - Thinning: None
    - K = 25
    - L = 5
    - `W*` ~ Gamma(1, rate=1/2)
    - `p_i` ~ Beta(1, 9) => `omega_q` ~ Normal(-3, 1.3)
    - phi = (0, 1, 10, 100) (command line arg)
    - Temperatures = (1, 1.003, 1.006, 1.01)
    - Save all states
    - Patients Data: (001_d31, 007_d35, 010_d35)
    - removed markers:
        - [2, 3, 4, 6, 8, 9, 10, 11, 12, 18, 19, 20, 21, 23, 26, 30, 31]
        - These markers are either mostly non-expressed or mostly expressed
    - Run with samples 1 and 2 only.
    - ~~Z: SB-IBP~~
    - Repulsive function = `pow(d, phi)`

- Visualizations
    - [ ] TSNE
        - Color by cluster

- Comparisons
    - Not in paper
        - [ ] Flowsom
        - [ ] Mclust
