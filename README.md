# PGE383-Thermodynamics
Repository for Fall 2022 course assignments in PGE383: Advanced Thermodynamics and Phase Behavior

To run the code: 
1. Change input_file.yml to fit problem specifications.
2. One needs only to run the main file for individual assignments. An output directory will be created with numerical outputs and any figures.

---
**For Tesleem:**

Work for Problems 1, 2a, 2b, and 3 are written in the attached HW4.pdf file. The code for 2c is contained in this repository.

I use YAML format for input files. PyYAML is a dependency for the code to run. You should be able to simply change the P and T in the input file and rerun HW4_main.py.

The code automatically creates an output directory if one does not already exist. Within the directory will be 3 types of files:
1) A .txt file with stdout outputs. This will contain the cubic equation coefficients, solutions to Cardano's method and the calculated molar volume.
2) A .png file with a plot of the cubic EoS function vs. Z for a graphical method of validating the code
3) An .xlsx file with the data from the cubic EoS function vs. Z. I created arrays in Python instead of directly using Excel but graphically found the root.

Equation of State and EoS utility functions are found in eos.py. Functions related to input/output are found in io_utils.py. Cardano's method function is found in solve.py.


