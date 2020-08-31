# eDITH-SamplingStrategy

MATLAB code and data supporting Carraro L, Stauffer JB, Altermatt F. How to design optimal eDNA sampling strategies for biomonitoring in river networks, *environmental DNA* (accepted), 2020.

## Main scripts

- `prepare_OCN.R`: creates and saves the OCN (`OCN_eDNA.rda`), and produces Figure 4 of the manuscript
- `perform_simulations.R`: reads `OCN_eDNA.rda` and runs simulations. Results are saved in `results_all.rda`
- `analyze_results.R`: reads `results_all.rda` and produces Figures 1, 3, 5, 6, S1, S2, S3, S4, S5, S6, S8, S9, S10, S11
- `prepare_OCN_SI.R`: creates and saves the OCN (`OCN_eDNA_SI.rda`) used in the supporting information analysis (with Q0 = 2.5 m<sup>3</sup>s<sup>-1</sup>).
- `perform_simulations_SI.R`: reads `OCN_eDNA_SI.rda` and runs simulations. Results are saved in `results_tau4_Q25.rda` (as an example, all other 'results' files can be produced as well).
- `analyze_results_SI.R`: reads `results_all.rda`, `results_tau4_Q25.rda`, `results_tau1_Q25.rda`, `results_tau1_Q10.rda`, and produces Figure S7.
- `eDITH_functions.R`: functions implementing the eDITH model