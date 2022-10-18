# PGE384-Thermodynamics
Repository for Fall 2022 course assignments in PGE384: Advanced Thermodynamics and Phase Behavior

To run the code: 
1. Change input_file.yml to fit problem specifications.
2. One needs only to run the main file for individual assignments. An output directory will be created with numerical outputs and any figures.

---
**For Tesleem:**

Work for Problems 1a, 1c, 2, and 3 are written in the attached HW5.pdf file. The code for 1c is contained in this repository.

I use YAML format for input files. PyYAML is a dependency for the code to run. You should be able to simply change the P and T in the input file and rerun HW5_main.py. I also specify my desired critical temperature, critical pressure, and acentric factor. I supply input files from the Input_Files/ directory. The input file for this assignment is titled "hw5_input.yml"

The code automatically creates an output directory if one does not already exist. For this assignment, it will be called HW5_Output/. Within the directory will be 3 types of files:
1) A .txt file with stdout outputs. This will contain the cubic equation coefficients, solutions to Cardano's method, the calculated molar volume, and the entropy and enthalpy departures.
2) A .png file with a plot of the cubic EoS function vs. Z for a graphical method of validating the code
3) An .xlsx file with the data from the cubic EoS function vs. Z. I created arrays in Python instead of directly using Excel but graphically found the root.

Equation of State and EoS utility functions are found in eos.py. Functions related to input/output are found in io_utils.py. Cardano's method function is found in solve.py. Departure functions for this assignment are in the pr_utils.py file. 


