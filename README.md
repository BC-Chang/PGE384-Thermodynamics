# PGE384-Thermodynamics
Repository for Fall 2022 course assignments in PGE384: Advanced Thermodynamics and Phase Behavior

To run the code: 
1. Change input_file.yml to fit problem specifications.
2. One needs only to run the main file for individual assignments. An output directory will be created with numerical outputs and any figures.

---
**For Tesleem:**

*Homework 8:*

The work for problem 1 is shown in the attached PDF file. 

Input Files: 
I use YAML format for input files. PyYAML is a dependency for the code to run. I use separate input files for problems 1 and 2, they can be found in the Input_Files/ directory. Setting ```Pvap = None``` will automatically initialize vapor pressure to the vapor pressure given by Wilson's correlation. You should be able to simply change the P and T in the input file and rerun HW8_main.py. I also specify the desired critical temperature, critical pressure, and acentric factor for each fluid as a list (state variables grouped together).

```HW8_main.py``` is split into 2 functions that can run independently. ```p2_main()``` corresponds to Problem 2 and ```p3_main()``` corresponds to Problem 3. One can independently run each problem by commenting/uncommenting the function calls within the script.

The code automatically creates an output directory if one does not already exist. For this assignment, it will be called HW8_Output/. There will be 2 subdirectories, Problem 2/ and Problem 3/. 

Within the HW8_Output/Problem 2/ directory will be 2 files:
1) A .txt file with stdout outputs. This will contain the final vapor pressure for each component, phase compositions, K values, Dew point and bubble point, and K values at the specified pressure and temperature.
2) A .png file titled *Problem2C_RR.png* with the Rachford-Rice vs. $\beta_v$ plot.

Within the HW8_Output/Problem 3/ directory will be 1 file:
1) A .txt file with stdout outputs. This will contain the fugacity coefficients and fugacity (in psi) of each component. 

Equation of State and EoS utility functions are found in ```eos.py```. Functions related to input/output are found in ```io_utils.py```. Cardano's method function and flash calculations are found in ```solve.py```. Departure functions and fugacity calculations for this assignment are in the ```pr_utils.py``` file. 

Changes since project 1:
 - Created a ```singlecomponent_utils.py``` function file that contains functions for single component flash and vapor pressure calculations
 - Created a ```multicomponent_utils.py``` function file that contains utilities for calculating bubble point, dew point, and phase compositions
 - Created a ```multicomponent_solve.py``` function file that contains a general-use function for performing Newton-Raphson iteration for root finding, a lambda function for the Rachford-Rice equation, and a function for finding the root of the Rachford-Rice equation.
 - Added a ```fugacity_coefficient_multicomponent()``` function for calculating the fugacity coefficient and fugacity (in Pa)
 - Created a ```Unit_Converter()``` class to help with unit conversions.




