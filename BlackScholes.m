%Define variables of problem
T=3; %final time
r=0.5; %interest rate of bank account
sigma=0.05; %volatility of stock
K=100; %strike price
S0=85; %initial stock price

%Calculate true option price
d1=(log(S0/K)+(r+0.5*sigma^2)*T)/(sigma*sqrt(T));
d2=d1-sigma*sqrt(T);
V0True=S0*normcdf(d1)-K*exp(-r*T)*normcdf(d2);
disp(['True value of call option is: ',num2str(V0True)]);

%Now use MC to estimate
%% True solution
h=[0.25,0.2,0.1,0.05,0.025];
MCError=zeros(1,length(h));
M=1e6; %number of samples to take
for j=1:length(h)
V0=0; %initialise
for i=1:M
    [~,S_t,~]=explicitEulerGBM(S0,r,sigma,T/h(j),T/h(j),0,T);
    payoff=S_t(length(S_t))-K; %payoff of call option
    if payoff>0
        V0=V0+(1/M)*payoff;
    end
end
V0=exp(-r*T)*V0;
MCError(j)=abs(V0True-V0);
end
varNames={'h','ErrorInOptionPriceEstimate'};
ErrorTable=table(h',MCError','VariableNames',varNames);
disp(ErrorTable);
%% Euler method 
h=[0.25,0.2,0.1,0.05,0.025];
MCErrorEuler=zeros(1,length(h));
M=1e6; %number of samples to take
tic
for j=1:length(h)
V0Euler=0; %initialise
for i=1:M
    [S_t,~,~]=explicitEulerGBM(S0,r,sigma,T/h(j),T/h(j),0,T);
    payoff=S_t(length(S_t))-K; %payoff of call option
    if payoff>0
        V0Euler=V0Euler+(1/M)*payoff;
    end
end
V0Euler=exp(-r*T)*V0Euler;
MCErrorEuler(j)=abs(V0True-V0Euler);
end
toc
varNames={'h','ErrorInOptionPriceEstimate'};
ErrorTable=table(h',MCErrorEuler','VariableNames',varNames);
disp(ErrorTable);
%% Milstein method
h=[0.25,0.2,0.1,0.05,0.025];
MCErrorMilstein=zeros(1,length(h));
M=1e6; %number of samples to take
tic
for j=1:length(h)
V0Milstein=0; %initialise
for i=1:M
    [S_t,~,~]=milsteinGBM(S0,r,sigma,T/h(j),T/h(j),0,T);
    payoff=S_t(length(S_t))-K; %payoff of call option
    if payoff>0
        V0Milstein=V0Milstein+(1/M)*payoff;
    end
end
V0Milstein=exp(-r*T)*V0Milstein;
MCErrorMilstein(j)=abs(V0True-V0Milstein);
end
toc
varNames={'h','ErrorInOptionPriceEstimate'};
ErrorTable2=table(h',MCErrorMilstein','VariableNames',varNames);
disp(ErrorTable2);
%% Weak Euler method
h=[0.25,0.2,0.1,0.05,0.025];
MCErrorWeakEuler=zeros(1,length(h));
M=1e6; %number of samples to take
tic
for j=1:length(h)
V0WeakEuler=0; %initialise
for i=1:M
    S_t=randomWalkGBM(S0,r,sigma,T/h(j),0,T);
    payoff=S_t(length(S_t))-K; %payoff of call option
    if payoff>0
        V0WeakEuler=V0WeakEuler+(1/M)*payoff;
    end
end
V0WeakEuler=exp(-r*T)*V0WeakEuler;
MCErrorWeakEuler(j)=abs(V0True-V0WeakEuler);
end
toc
varNames={'h','ErrorInOptionPriceEstimate'};
ErrorTable3=table(h',MCErrorWeakEuler','VariableNames',varNames);
disp(ErrorTable3);

%% Talay-Tubaro method of order two
h=[1,0.5,0.25,0.125,0.0625];
MCErrorTalay=zeros(1,length(h));
M=1e7; %number of samples to take
tic
for j=1:length(h)
V0Talay=0; %initialise
for i=1:M
    [extrapValue,~,~]=TalayTubaroOrderTwo(S0,r,sigma,T/h(j),0,T);
    payoff=extrapValue-K; %payoff of call option
    if payoff>0
        V0Talay=V0Talay+(1/M)*payoff;
    end
end
V0Talay=exp(-r*T)*V0Talay;
MCErrorTalay(j)=abs(V0True-V0Talay);
end
toc
varNames={'h','ErrorInOptionPriceEstimate'};
ErrorTableTalay=table(h',MCErrorTalay','VariableNames',varNames);
disp(ErrorTableTalay);
%% Talay-Tubaro method of order three 
h=[3,1.5,0.75,0.375,0.1875];
MCErrorTalay2=zeros(1,length(h));
M=1e6; %number of samples to take
tic
for j=1:length(h)
V0Talay2=0; %initialise
for i=1:M
    extrapValue2=TalayTubaroOrderThree(S0,r,sigma,T/h(j),0,T);
    payoff=extrapValue2-K; %payoff of call option
    if payoff>0
        V0Talay2=V0Talay2+(1/M)*payoff;
    end
end
V0Talay2=exp(-r*T)*V0Talay2;
MCErrorTalay2(j)=abs(V0True-V0Talay2);
end
toc
varNames={'h','ErrorInOptionPriceEstimate'};
ErrorTableTalay2=table(h',MCErrorTalay2','VariableNames',varNames);
disp(ErrorTableTalay2);
