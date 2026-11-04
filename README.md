# distributed-i-IOSS-verification
MATLAB implementation for distributed verification of exponential incremental input/output-to-state stability for a nonlinear train model.


This repository contains the code of the (not yet published) paper

Christian Gatke, Julian D. Schiller, Matthias A. Müller. Small-gain analysis of incremental input/output-to-state stability for large-scale distributed systems. 2026.


The following MATLAB files are provided:
- `distributed_verification/D1_train_simulation.m`: Simulation of a train acceleration process.
- `distributed_verification/D2_1_LMI_SDP_outer.m`: LMI verification for the boundary carriages.
- `distributed_verification/D2_2_LMI_SDP_inner.m`: LMI verification for the intermediate carriages.
- `distributed_verification/D3_1_LMI_check_outer.m`: Check satisfaction of LMI conditions for the boundary carriages.
- `distributed_verification/D3_2_LMI_check_inner.m`: Check satisfaction of LMI conditions for the intermediate carriages.
- `distributed_verification/D4_SmallGain_condition.m`: Setting up gain matrices and evaluation of their spectral radii.
- `Centralized_verification.m`: Centralized verification of the train model for different numbers of carriages.

Distributed verification:
1. Run `D1_train_simulation.m` to simulate the acceleration process and determine a suitable upper bound on the relative velocities between carriages.
2. Select this bound in `D2_1_LMI_SDP_outer.m` and `D2_2_LMI_SDP_inner.m` and perform the LMI verification.
3. Load the resulting .mat-files into `D3_1_LMI_check_outer.m` and `D3_2_LMI_check_inner.m` to check whether the LMI conditions are indeed satisfied.
4. If the LMI conditions are satisfied, load the corresponding .mat-files in `D4_SmallGain_condition.m` to verify the small-gain condition.

Centralized verification:
Run `Centralized_verification.m` to perform the centralized verification for a specified range of carriage numbers and analyze the computation time.


The results described in the paper were achieved using the follwoing software/hardware configurations:
- MATLAB, Version 24.2.0.2923080 (R2024b) Update 6
- YALMIP, Version 20250626
- MOSEK, Version 11.0.28
- Windows 11, 12th Gen Intel(R) Core(TM) i7-1270P, 16 GB RAM