# StochasticNumericsDissertation
The MATLAB code for fourth-year dissertation project on numerical methods for stochastic differential equations.
These codes were used to produce the graphs and tables in the dissertation report. 
This file contains a description of which codes are relevant to which chapter, and also informs the reader how to use them. 
Open the files for a description of what they do and how they work.

Chapter 3 (Brownian Motion - Simulation in 1D): 
The script BrownianSimulation.m is used to produce figure 4 (page 41). 
This code is self-contained and does not make calls to any other functions, thus can be ran as is by the reader.

Chapter 5 (Stochastic Numerics: Strong Taylor Approximations - Implementing the Explicit Euler scheme):
The function explicitEulerGBM is called from code sections 1,2 and 3 within GBMDriver.m to produce figure 6/7 (page 72/73).
Use the 'Run Section' feature to run this code.

Chapter 5 (Stochastic Numerics: Strong Taylor Approximations - Generalising the implementation of the Explicit Euler Method):
The function explicitEulerMethod.m, which calls the function evaluateFunctions.m, is called from code sections 4,5,6,7 and 9 
(evaluateFunctions.m is called from explicitEulerMethod.m through GBMDriver.m in sections 4,6 and 7) within GBMDriver.m to 
produce figure 8/9 (page 76/77). Note that in sections 4 and 5, there is a function test2 used as an input of explicitEulerMethod.m.
To run these sections, the reader will need to copy and paste each section in turn to the bottom of the script, otherwise
an error will be displayed (functions must be written at the bottom of script files).
Use the 'Run Section' feature to run this code. 
