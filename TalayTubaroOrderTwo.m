function extrapOrderTwo = TalayTubaroOrderTwo(initialS,mu,sigma,N_approx,t0,T)
%TalayTubaroOrderTwo.m: Function taking inputs defined for the problem
%along with a number of time-steps, which calls randomWalkGBM.m twice to
%obtain estimates before outputting the improved estimate with second-order
%accuracy.

u_h1=randomWalkGBM(initialS,mu,sigma,N_approx,t0,T); %get first approximation by weak Euler method 
u_h2=randomWalkGBM(initialS,mu,sigma,N_approx*2,t0,T); %get second approximation with twice number of time steps as above
h1_approx=u_h1(length(u_h1)); %get the approximation of process at final time
h2_approx=u_h2(length(u_h2)); %get the approximation of process at final time
h1=(T-t0)/N_approx; %define the two time steps
h2=0.5*h1;
C0=(h1_approx-h2_approx)/(h2-h1);
extrapOrderTwo=h1_approx+C0*h1; %form the improved estimate 



