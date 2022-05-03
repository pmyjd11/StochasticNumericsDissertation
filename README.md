# StochasticNumericsDissertation
The MATLAB code for fourth-year dissertation project on numerical methods for stochastic differential equations.
These codes were used to produce the graphs and tables in the dissertation report. 
This file contains a description of which codes are relevant to which chapter, and also informs the reader how to use them. 
Open the files for a description of what they do and how they work.

Chapter 3 (Brownian Motion - Simulation in 1D): 
The script BrownianSimulation.m is used to produce figure 4 (page 41). 
This code is self-contained and does not make calls to any other functions, thus can be ran as is by the reader.

Chapter 5 (Stochastic Numerics: Strong Taylor Approximations - Implementing the Explicit Euler scheme):
The function explicitEulerGBM is called from code section 1 within GBMDriver.m to produce figure 6/7 (page 72/73).
Use the 'Run Section' feature to run this code.

Chapter 5 (Stochastic Numerics: Strong Taylor Approximations - Generalising the implementation of the Explicit Euler Method):
The function explicitEulerMethod.m, which calls the function evaluateFunctions.m, is called from code sections 4 and 5  
(evaluateFunctions.m is called from explicitEulerMethod.m through GBMDriver.m in sections 4) within GBMDriver.m to 
produce figure 8/9 (page 76/77). Note that in sections 4 and 5, there is a function test2 used as an input of explicitEulerMethod.m.
To run these sections, the reader will need to copy and paste each section in turn to the bottom of the script, otherwise
an error will be displayed (functions must be written at the bottom of script files).
Use the 'Run Section' feature to run this code. 

Chapter 5 (Stochastic Numerics: Strong Taylor Approximations - Analysing the Euler scheme: The absolute error criterion):
The function explicitEulerGBM.m is called from code section 2 within GBMDriver.m to produce the table on page 82.
The function explicitEulerMethod.m, which calls the function evaluateFunctions.m, is called from code sections 6 and 7 
(evaluateFunctions.m is called from explicitEulerMethod.m through GBMDriver.m in section 6) within GBMDriver.m to 
produce figure 10/11 (page 84).
Use the 'Run Section' feature to run this code.

Chapter 5 (Stochastic Numerics: Strong Taylor Approximations - Milstein scheme):
The function milsteinGBM.m is called from code section 8 within GBMDriver.m to produce figure 13 on page 90, and it is 
called from code section 11 to produce the table on page 91. Figure 14 (page 92) is produced through running code section 
12 of GBMDriver.m.
Use the 'Run Section' feature to run this code. 


