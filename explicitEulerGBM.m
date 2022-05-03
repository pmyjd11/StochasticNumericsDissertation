function [eulApprox,trueGBM,samplePath] = explicitEulerGBM(initialS,mu,sigma,N_approx,N_true,t0,T,zValues)
%explicitEulerGBM: Function that takes inputs from the problem, two numbers
%of time-steps and also an optional vector zValues which holds the sample
%path (if already generated), and outputs the Euler approximation, true
%solution and the sample path that generates them. 

if nargin==7 %If values of sample path not given then generated here
    zValues=zeros(1,N_true); %to hold the sample path
    for i=1:N_true
        zValues(i)=normrnd(0,1); %generate sample path
    end
end
samplePath=zValues; %save this as an output so that approximations on same sample path can be obtained

%First get a trajectory of the true solution
trueGBM=zeros(1,N_true+1); 
trueGBM(1)=initialS;
hTrue=(T-t0)/N_true; %find size of time-step for true 
t=t0; %this will keep track of the time for true solution
for i=1:N_true
    t=t+hTrue; %update this for the true solution
    trueGBM(i+1)=trueGBM(i)*exp((mu-0.5*sigma^2)*hTrue+sigma*sqrt(hTrue)*zValues(i)); %calculate next element in true solution
end

if N_true<N_approx
    disp('Error: You must input a larger number of time-steps for the true solution than the approximation');  
elseif mod(N_true,N_approx)~=0 %ensures we can use the same sample path properly
    N_approx=N_true/floor(N_true/N_approx);
    disp('Warning: Number of time-steps for approximation not a factor of number of time-steps for exact');
    disp(['Number of time-steps in approximation changed to N=',num2str(N_approx)]);
end


hApprox=(T-t0)/N_approx; %find the size of each time-step for approximate 
eulApprox=zeros(1,N_approx+1); %initialise the output vectors
eulApprox(1)=initialS; %assign the first value in each vector to initialS
factor=hApprox/hTrue; %used to match the sample paths of exact with approximate

for i=1:N_approx
    increment=0; %reset value
    for k=1:factor
        increment=increment+sqrt(hTrue)*zValues(k+(i-1)*factor); %sum up to obtain W(t+hApprox)-W(t) for discretisation points t
    end
    eulApprox(i+1)=eulApprox(i)+hApprox*mu*eulApprox(i)+sigma*eulApprox(i)*increment; %calculate next element in approximation
end



