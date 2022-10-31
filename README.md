# PGE384-Thermodynamics
Repository for Fall 2022 course assignments in PGE384: Advanced Thermodynamics and Phase Behavior

To run the code: 
1. Change input_file.yml to fit problem specifications.
2. One needs only to run the main file for individual assignments. An output directory will be created with numerical outputs and any figures.

---
**For Tesleem:**

*Project 1:*

The work for problem 1a and 2a are shown in the attached PDF file. 

Input Files: 
I use YAML format for input files. PyYAML is a dependency for the code to run. I use separate input files for problems 1 and 2, they can be found in the Input_Files/ directory. Setting ```P = None``` will automatically initialize pressure to the vapor pressure given by Wilson's correlation. You should be able to simply change the P and T in the input file and rerun Project1_main.py. I also specify the desired critical temperature, critical pressure, and acentric factor for the working fluid.

The code automatically creates an output directory if one does not already exist. For this assignment, it will be called Project1_Output/. There will be 2 subdirectories, Problem 1/ and Problem 2/. 

Within the Project1_Output/Problem 1/ directory will be 5 files:
1) A .txt file with stdout outputs. This will contain the final vapor pressure and the equilibrium molar volumes at the specified temperature.
2) A .png file titled *PV_isotherm_<temp>.png* with a plot of PV isotherm at the specified temperature.
3) Two .xlsx files with the pressure and molar volume values from the PT flash calculation, used to plot the combined PV isotherm. Additional .xlsx files will be generated if a temperature other than 313.15K (40C) and 343.15 (70C) is used.
4) A .png file titled "Combined_PV_isotherm" with the PV isotherm values in each .xlsx file, along with the critical point.

Within the Project1_Output/Problem 2/ directory will be 4 files:
1) A .txt file with stdout outputs. This will contain the vapor pressure, equilibrium molar volumes at specified temperature, as well as answers to the questions posed in problems 2c and 2d. 
2) A .png file with the PV Isotherm and Molar Gibbs Free Energy curves.
3) A .png file stable, metastable, and unstable regions in the Molar Gibbs Free Energy curves.
4) A .png file with the stable regions of the Molar Gibbs Free Energy curve.

Equation of State and EoS utility functions are found in ```eos.py```. Functions related to input/output are found in ```io_utils.py```. Cardano's method function and flash calculations are found in ```solve.py```. Departure functions and fugacity calculations for this assignment are in the ```pr_utils.py``` file. 


