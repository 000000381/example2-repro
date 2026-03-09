# example2-repro
Reproducibility code for Example 2 in the revised manuscript.

Tested environment: MATLAB R2024a (not tested on other versions).

How to run:
1.Download Run_Example.m.
2.Open MATLAB R2024a and set the current folder to the repository root.
3.Run Run_Example2.m.

What it does:
Simulates Example 2 in the revised manuscript, and reproduces the main plots reported in the paper.

Notes:
The file is self-contained: helper functions (dynamic3, dM_derive) are included at the end of the same .m file.
Parameters are written explicitly for readability; they can be moved outside the post-processing loop for efficiency without changing the simulation results.
