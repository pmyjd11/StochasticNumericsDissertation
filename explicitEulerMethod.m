function [eulApprox,samplePath] = explicitEulerMethod(funcIn,params,x0,N,t0,T,zValues,factor)
%explicitEulerMethod.m: Function taking inputs of a function funcIn (used
%to evaluate the coefficients of an SDE), params (parameters which may be
%an input of funcIn), the initial state of the process, a number of
%time-steps N, the interval limits, optional sample path vector zValues and
%an optional integer factor which tells us how many times N goes into any 
%N used in a previous call to the function. The outputs include the Euler
%approximation and the sample path used to generate this approximation.

h=(T-t0)/N; %define step-size
eulApprox=zeros(1,N+1); %initialise vector to fill with the approximations
eulApprox(1)=x0; %set the first element equal to the initial state x0
t=t0; %keep track of the time
evals=funcIn(0,0,params,0); %this is so we can obtain the dimension of the problem
if nargin==6 % if zValues and factor not passed to function
    factor=1; %set so that zValues generated for N time-steps
    zValues=zeros(1,(length(evals)-1)*N);  %initialise
    for k=1:length(zValues)
        zValues(k)=normrnd(0,1); %generate samples from standard normal distribution
    end
end
hTrue=(T-t0)/(factor*N); %this only important if sample path passed to function in zValues
samplePath=zValues;
evals=funcIn(t,eulApprox(1),params,zValues); %call for function evaluations, also tells us the dimension of the problem
for i=1:N
    eulApprox(i+1)=eulApprox(i)+h*evals(1); %add previous and dt 
    for r=2:length(evals) %for loop adds the brownian part
        increment=0; %reset value
        for k=1:factor %add up the increments 
        increment=increment+sqrt(hTrue)*zValues(r-1+(length(evals)-1)*(k-1)+(i-1)*(length(evals)-1)*factor); %sum up to obtain W(t+hApprox)-W(t) for discretisation points t
        end
        eulApprox(i+1)=eulApprox(i+1)+evals(r)*increment;
    end
    t=t+h; %update the time
    evals=funcIn(t,eulApprox(i+1),params,zValues); %call for function evaluations
end


