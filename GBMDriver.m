%% Produce plots showing true trajectory with approximations for various h
%Code calls the explicitEulerGBM function to compute the true GBM solution
%along with an approximation computed using the explicit Euler method. It
%then produces a plot of the two along with the solution of the
%deterministic equation (i.e. where sigma=0)
mu=0.08; %declare mu variable
N=8; %initialise number of time-steps for approximation
[eulApprox,trueGBM,samplePath]=explicitEulerGBM(160,mu,0.08,N,2048,0,8); %call function to compute true solution and approximation and return sample path
tApprox=linspace(0,8,N+1); %make a vector to allow us to plot approximate solution
t=linspace(0,8,2049);
plot(t,trueGBM,'b-'); %plot true solution
hold on
plot(tApprox,eulApprox,'r-'); %plot approximate true solution
for i=1:3
    N=2*N; %double number of time-steps in approximation
    tApprox=linspace(0,8,N+1); 
    [eulApprox,~,~]=explicitEulerGBM(160,mu,0.08,N,2048,0,8,samplePath); %call function, providing same sample path as above
    plot(tApprox,eulApprox);
end
hold off 
xlabel('t');
ylabel('S(t)');
legend('True (N=2048)','Euler (N=8)','Euler (N=16)','Euler (N=32)','Euler (N=64)','Location','northwest');

%% Analyse the absolute error criterion
%This section of code takes 100 samples of various numbers of time-steps N
%and estimates the absolute error at T=8 using the arithmetic mean, which
%is an unbiased estimator. 

mu=0.08; %declare mu variable
N=4; %initialise number of time-steps for approximation
absoluteErrorMean=zeros(1,10); %initialise
stepsize=zeros(1,10);
sampleNumber=500; %number of samples to take for each N
for j=1:10
    stepsize(j)=8/N;
    absoluteError=zeros(1,sampleNumber); %reset this vector
    for i=1:sampleNumber
        [eulApprox,trueGBM]=explicitEulerGBM(160,mu,0.08,N,2048,0,8); %call function to compute true solution and approximation
        absoluteError(i)=abs(trueGBM(length(trueGBM))-eulApprox(length(eulApprox)));
    end
    absoluteErrorMean(j)=mean(absoluteError); 
    N=2*N;
end
varNames=["h","AbsoluteErrorAtT8"];
absoluteErrorTable=table(stepsize',absoluteErrorMean','VariableNames',varNames);
disp(absoluteErrorTable);
root2ratio=zeros(1,10); %now look at ratios
root2ratio(1)=absoluteErrorMean(1);
trueRatio=zeros(1,10);
trueRatio(1)=1;
for i=1:9
    root2ratio(i+1)=root2ratio(i)/sqrt(2); %if the ratio between the absolute errors is sqrt(2) these are the absolute errors we would see
    trueRatio(i+1)=absoluteErrorMean(i)/absoluteErrorMean(i+1); %this is the true ratio between the absolute errors given by the estimator
end
varNames=["h","AbsoluteErrorAtT8","AbsoluteErrorUnderRoot2Ratio","TrueRatio"];
errorTableWithRatio=table(stepsize',absoluteErrorMean',root2ratio',trueRatio','VariableNames',varNames);
disp(errorTableWithRatio);

%% Confidence intervals for the absolute error criterion

Napprox=64; %number of time-steps in the approximation
M=[10,20,40,100]; %number of batches
N=100; %number of simulations in each batch
meanEps=zeros(1,length(M));
varEps=zeros(1,length(M));
standardErrorEps=zeros(1,length(M));
confidenceIntervalEps=zeros(length(M),2);
for i=1:length(M)
    epsJ=zeros(1,M(i));
for j=1:M(i)
    for k=1:N
        [eulApprox,trueGBM]=explicitEulerGBM(160,mu,0.08,Napprox,2048,0,8); %call function to compute true solution and approximation
        epsJ(j)=epsJ(j)+(1/N)*abs(trueGBM(length(trueGBM))-eulApprox(length(eulApprox)));
    end
end
meanEps(i)=mean(epsJ);
varEps(i)=0;
for j=1:M(i)
    varEps(i)=varEps(i)+(1/(M(i)-1))*(epsJ(j)-meanEps(i))^2;
end
standardErrorEps(i)=tinv(0.95,M(i)-1)*sqrt(varEps(i)/M(i));
confidenceIntervalEps(i,1)=meanEps(i)-standardErrorEps(i);
confidenceIntervalEps(i,2)=meanEps(i)+standardErrorEps(i);
end
errorbar(M,meanEps,standardErrorEps,'s','MarkerEdgeColor','red','MarkerFaceColor','red');
xlim([0,110]);
ylim([0.9*confidenceIntervalEps(1,1),1.1*confidenceIntervalEps(1,2)]);
xlabel('Number of batches (M)');
ylabel('Absolute error');
title('Confidence intervals with h=7.8125e-3');

%% Use the more general function for the approximation
%This section of code tests the more general function explicitEulerMethod.m
%and the function evaluateFunctions.m which is called through the former.
%First it is applied to GBM then on the SDE df=Wdt+tdW which has known
%solution.
[~,trueGBM,samplePath]=explicitEulerGBM(160,0.08,0.08,32,1024,0,8); %call function to compute true solution
[eulApprox,~]=explicitEulerMethod(@evaluateFunctions,[0.08,0.08],160,32,0,8,samplePath,32); %call function to compute approximation
t=linspace(0,8,1025); %declare t vector so we can plot true solution
tApprox=linspace(0,8,33); %declare so we can plot the approximation
subplot(2,1,1); %rest of this part plots GBM solution/approximation
plot(t,trueGBM,'b-');
hold on
plot(tApprox,eulApprox,'r-'); 
hold off 
xlabel('t');
ylabel('S(t)');
legend('True (sigma=0.08,N=1024)','Euler approx. (sigma=0.08,N=32)','Location','northwest');
title('Geometric Brownian Motion');

trueSol=zeros(1,1025); 
trueSol(1)=0; %initial state 
h=8/1024;
time=0;
tApprox=linspace(0,8,33);
zValues=normrnd(0,1,[1,1024]);
for i=1:1024 %for loop to find the true solution
    time=time+h; %update time
    increment=0;
    for j=1:i
        increment=increment+sqrt(h)*zValues(j);
    end
    trueSol(i+1)=time*increment;
end
[eulApprox,~]=explicitEulerMethod(@test2,[8/32,0],0,32,0,8,zValues,32); %call again to find approximation to second SDE
subplot(2,1,2); %rest of this section plots the above
plot(t,trueSol,'b-');
hold on
plot(tApprox,eulApprox,'r-');
hold off
xlabel('t');
ylabel('f(t)');
legend('True (N=1024)','Euler approx. (N=32)','Location','best');
title('f(t,W(t))=tW(t)');

%function evals=test2(t,~,params,zValues) %function to evaluate coefficients in df=W(t)dt+tdW(t)
%evals=zeros(1,2);
%if zValues~=0
%for k=1:t/params(1)
%evals(1)=evals(1)+sqrt(params(1))*zValues(k);
%end
%evals(2)=t;
%end
%end


%% Repeat second problem for various values of h

trueSol=zeros(1,1025); 
trueSol(1)=0; %initial state 
h=8/1024;
time=0;
zValues=normrnd(0,1,[1,1024]);
for i=1:1024 %for loop to find the true solution
    time=time+h; %update time
    increment=0;
    for j=1:i
        increment=increment+sqrt(h)*zValues(j);
    end
    trueSol(i+1)=time*increment;
end
plot(t,trueSol,'b-');
hold on
N=128;
for i=1:3
    tApprox=linspace(0,8,N+1);
    zValues(1025)=8/N; %add in the value of h for each
    [eulApprox,~]=explicitEulerMethod(@test2,[8/N,0],0,N,0,8,zValues,1024/N); %call again to find approximation to second SDE
    plot(tApprox,eulApprox);
    N=2*N;
end
hold off 
xlabel('t');
ylabel('f(t)');
legend('True (N=1024)','Euler approx. (N=128)','Euler approx. (N=256)','Euler approx. (N=512)','Location','best');
title('f(t,W(t))=tW(t)');

%function evals=test2(t,~,params,zValues) %function to evaluate coefficients in df=W(t)dt+tdW(t)
%evals=zeros(1,2);
%if zValues~=0
%for k=1:t/params(1)
%evals(1)=evals(1)+sqrt(params(1))*zValues(k);
%end
%evals(2)=t;
%end
%end

%% Mean error consideration
Napprox=64; %number of time-steps in the approximation
M=[10,20,40,100]; %number of batches
N=100; %number of simulations in each batch
meanEps=zeros(1,length(M));
varEps=zeros(1,length(M));
standardErrorEps=zeros(1,length(M));
confidenceIntervalEps=zeros(length(M),2);
averageXT=160*exp(0.08*8); %this is the average value of the GBM problem
for i=1:length(M)
    epsJ=zeros(1,M(i));
for j=1:M(i)
    for k=1:N
        [eulApprox,~]=explicitEulerMethod(@evaluateFunctions,[0.08,0.08],160,Napprox,0,8); %call function to give us an approximation
        epsJ(j)=epsJ(j)+(1/N)*eulApprox(length(eulApprox)); %form arithmetic mean of trajectories at final time T
    end
    epsJ(j)=epsJ(j)-averageXT; %subtract average of true solution for mean error of batch
end
meanEps(i)=mean(epsJ); %find average
varEps(i)=0;
for j=1:M(i)
    varEps(i)=varEps(i)+(1/(M(i)-1))*(epsJ(j)-meanEps(i))^2; %this for loop finds the variance
end
standardErrorEps(i)=tinv(0.95,M(i)-1)*sqrt(varEps(i)/M(i));
confidenceIntervalEps(i,1)=meanEps(i)-standardErrorEps(i);
confidenceIntervalEps(i,2)=meanEps(i)+standardErrorEps(i);
end
figure(1);
errorbar(M,meanEps,standardErrorEps,'s','MarkerEdgeColor','red','MarkerFaceColor','red');
xlim([0,110]);
ylim([-10,10]);
xlabel('Number of batches (M)');
ylabel('Mean error');
title('Confidence intervals with h=7.8125e-3');

%% Repeat but vary step-size
NForApproximation=4; %number of time-steps in the approximation
M=100; %number of batches
N=100; %number of simulations in each batch
meanEps=zeros(1,10);
varEps=zeros(1,10);
standardErrorEps=zeros(1,10);
confidenceIntervalEps=zeros(10,2);
averageXT=160*exp(0.08*8); %this is the average value of the GBM problem
stepsize=zeros(1,10);
for i=1:10
    epsJ=zeros(1,100);
    stepsize(i)=8/NForApproximation;
for j=1:M
    for k=1:N
        [eulApprox,~]=explicitEulerMethod(@evaluateFunctions,[0.08,0.08],160,NForApproximation,0,8); %call function to give us an approximation
        epsJ(j)=epsJ(j)+(1/N)*eulApprox(length(eulApprox)); %form arithmetic mean of trajectories at final time T
    end
    epsJ(j)=epsJ(j)-averageXT; %subtract average of true solution for mean error of batch
end
meanEps(i)=mean(epsJ); %find average
varEps(i)=0;
for j=1:100
    varEps(i)=varEps(i)+(1/99)*(epsJ(j)-meanEps(i))^2; %this for loop finds the variance
end
standardErrorEps(i)=tinv(0.95,99)*sqrt(varEps(i)/100);
confidenceIntervalEps(i,1)=meanEps(i)-standardErrorEps(i);
confidenceIntervalEps(i,2)=meanEps(i)+standardErrorEps(i);
NForApproximation=2*NForApproximation;
end

figure(1);
errorbar(stepsize([1:4,10]),meanEps([1:4,10]),standardErrorEps([1:4,10]),'s','MarkerEdgeColor','red','MarkerFaceColor','red');
xlim([0,2]);
ylim([-20,20]);
xlabel('Step-size');
ylabel('Mean error');
title('Confidence intervals for varying h, M=N=100');

figure(2);
loglog(stepsize([1:4,10]),meanEps([1:4,10]),'-*','MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
xref = 10.^(-1:.1:0.5);
yref = (-1).*xref;
plot(xref,yref,'--b');
hold off
xlabel('h');
ylabel('Mean error');
title('Log-log plot of mean error against step-size');
xlim([0,2]);
ylim([-15,0]);
legend('Mean error','Reference line of slope 1','Location','best');

%% Milstein approximations for GBM problem
mu=0.08; %declare mu variable
N=8; %initialise number of time-steps for approximation
[milsteinApprox,trueGBM,samplePath]=milsteinGBM(160,mu,0.08,N,2048,0,8); %call function to compute true solution and approximation and return sample path
tApprox=linspace(0,8,N+1); %make a vector to allow us to plot approximate solution
t=linspace(0,8,2049);
plot(t,trueGBM,'b-'); %plot true solution
hold on
plot(tApprox,milsteinApprox,'r-'); %plot approximate true solution
for i=1:3
    N=2*N; %double number of time-steps in approximation
    tApprox=linspace(0,8,N+1); 
    [milsteinApprox,~,~]=milsteinGBM(160,mu,0.08,N,2048,0,8,samplePath); %call function, providing same sample path as above
    plot(tApprox,milsteinApprox);
end
hold off 
xlabel('t');
ylabel('S(t)');
legend('True (N=2048)','Milstein (N=8)','Milstein (N=16)','Milstein (N=32)','Milstein (N=64)','Location','northwest');

%% Comparison of Milstein with Euler
mu=0.08; %declare mu variable
N=128; %initialise number of time-steps for approximation
for i=1:2
    subplot(2,1,i);
[eulApprox,trueGBM,samplePath]=explicitEulerGBM(160,mu,0.08,N,2048,0,50); %call function to compute true solution and approximation and return sample path
[milsteinApprox,~,~]=milsteinGBM(160,mu,0.08,N,2048,0,50,samplePath); 
tApprox=linspace(0,50,N+1); %make a vector to allow us to plot approximate solution
t=linspace(0,50,2049);
plot(t,trueGBM,'b-'); %plot true solution
hold on
plot(tApprox,eulApprox,'r-');
plot(tApprox,milsteinApprox,'g-');
hold off
xlabel('t');
ylabel('S(t)');
legend('True (N=2048)','Euler (N=128)','Milstein (N=128)','Location','northwest');
end

%% Analysis of error in Milstein scheme for GBM (batch sizes)
Napprox=64; %number of time-steps in the approximation
M=[10,20,40,100]; %number of batches
N=100; %number of simulations in each batch
meanEps=zeros(1,length(M));
varEps=zeros(1,length(M));
standardErrorEps=zeros(1,length(M));
confidenceIntervalEps=zeros(length(M),2);
for i=1:length(M)
    epsJ=zeros(1,M(i));
for j=1:M(i)
    for k=1:N
        [milsteinApprox,trueGBM]=milsteinGBM(160,mu,0.08,Napprox,2048,0,8); %call function to compute true solution and approximation
        epsJ(j)=epsJ(j)+(1/N)*abs(trueGBM(length(trueGBM))-milsteinApprox(length(milsteinApprox)));
    end
end
meanEps(i)=mean(epsJ);
varEps(i)=0;
for j=1:M(i)
    varEps(i)=varEps(i)+(1/(M(i)-1))*(epsJ(j)-meanEps(i))^2;
end
standardErrorEps(i)=tinv(0.95,M(i)-1)*sqrt(varEps(i)/M(i));
confidenceIntervalEps(i,1)=meanEps(i)-standardErrorEps(i);
confidenceIntervalEps(i,2)=meanEps(i)+standardErrorEps(i);
end
errorbar(M,meanEps,standardErrorEps,'s','MarkerEdgeColor','red','MarkerFaceColor','red');
xlim([0,110]);
ylim([0.9*confidenceIntervalEps(1,1),1.1*confidenceIntervalEps(1,2)]);
xlabel('Number of batches (M)');
ylabel('Absolute error');
title('Confidence intervals with h=0.125');
%% Absolute error table with 500 simulations
Napprox=4; %number of time-steps in the approximation
N=500; %number of simulations in each batch
error=zeros(1,9);
stepsize=zeros(1,9);
for i=1:9
    for k=1:N
        [milsteinApprox,trueGBM]=milsteinGBM(160,mu,0.08,Napprox,2048,0,8); %call function to compute true solution and approximation
        error(i)=error(i)+(1/N)*abs(trueGBM(length(trueGBM))-milsteinApprox(length(milsteinApprox)));
    end
    stepsize(i)=8/Napprox;
    Napprox=2*Napprox;
end
ratio=zeros(1,9);
ratio(1)=1;
for i=1:8
    ratio(i+1)=error(i)/error(i+1);
end
MilsteinErrorTable=table(stepsize',error',ratio');
MilsteinErrorTable.Properties.VariableNames={'h','AbsoluteError','Ratio'};
disp(MilsteinErrorTable);
%% Plot of absolute error confidence intervals
Napprox=64; %number of time-steps in the approximation
N=100; %number of simulations in each batch
M=100;
epsJ=zeros(5,M);
meanEps=zeros(1,5);
varEps=zeros(1,5);
standardErrorEps=zeros(1,5);
confidenceIntervalEps=zeros(2,5);
stepsize=zeros(1,5);
for i=1:5 %fix number of time-steps
    for m=1:M %fix batch number
    for k=1:N %take N samples for each batch
        [milsteinApprox,trueGBM]=milsteinGBM(160,mu,0.08,Napprox,2048,0,8); %call function to compute true solution and approximation
        epsJ(i,m)=epsJ(i,m)+(1/N)*abs(trueGBM(length(trueGBM))-milsteinApprox(length(milsteinApprox))); %form average of the batch
    end        
    end
    meanEps(i)=mean(epsJ(i,1:M)); %take average of all the batches
    for j=1:M %find the variance
        varEps(i)=varEps(i)+(1/(M-1))*(epsJ(i,j)-meanEps(i))^2;
    end
    standardErrorEps(i)=tinv(0.95,M-1)*sqrt(varEps(i)/M);
    confidenceIntervalEps(1,i)=meanEps(i)-standardErrorEps(i);
    confidenceIntervalEps(2,i)=meanEps(i)+standardErrorEps(i);
    stepsize(i)=8/Napprox;
    Napprox=2*Napprox;
end
errorbar(stepsize,meanEps,standardErrorEps,'s','MarkerEdgeColor','red','MarkerFaceColor','red','MarkerSize',1.5);
xlabel('Step-size');
ylabel('Absolute error');
title('Confidence intervals for various step-size (M=N=100)');

