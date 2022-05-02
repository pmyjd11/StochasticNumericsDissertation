function evals = evaluateFunctions(t,x,params,zValues)
%evaluateFunctions.m: Function which takes inputs of t and x as well as a
%vector of further parameters params and uses them to find b(t,x) and
%sigma_r(t,x) within the SDE. The evaluations are placed in a single vector
%which can be used within the function that finds the approximate solution
%of the SDE.

evals=zeros(1,2); %initialise output vector, to be filled below
evals(1)=params(1)*x; %calculate mu*x
evals(2)=params(2)*x; %calculate sigma*x
%...
%Any further functions to evaluate parts of SDE go here


