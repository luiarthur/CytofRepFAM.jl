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

- Visualizations
    - [ ] TSNE
        - Color by cluster

- Comparisons
    - Not in paper
        - [ ] Flowsom
        - [ ] Mclust

- **This run is potentially good enough for the paper**
