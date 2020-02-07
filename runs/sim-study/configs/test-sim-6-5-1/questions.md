# Questions

## Bad:
- seed 3, scale 10 (why is the LPML bad?)
    - see Kmcmc 7 vs Kmcmc 8
    - mus and sig2 are different.
    - what does R.pdf look like?
- seed 1, scale 10 (why is the LPML bad?)
    - Notice the LPML looks bad. 
    - But this is explained by the graphs in R. For the
      Kmcmc values where there is a spike in LPML, the correct number of
      features were identified for exactly one sample. In no other Kmcmc
      runs were the number of features correctly identified.

## Good:
- seed 1, scale 0 (why is the LPML good?)
    - LPML plateaus at Kmcmc=7. But the correct Z is not identified.
      Two identical columns are present (and split evenly). This happened in
      Kmcmc=8 as well. So while the LPML looks good, the simulation truth was
      not really recovered.

## Maybe Good?:
- seed 1, scale 1 (why is the LPML bad?)
    - The LPML looks jagged. Difficult to choose Kmcmc.
    - notice that in R.pdf, when Kmcmc = 8, R1 = 6 and R2 = 5, which explains
      the high LPML.
