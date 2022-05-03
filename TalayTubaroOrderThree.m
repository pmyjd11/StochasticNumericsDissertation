function extrapOrderThree = TalayTubaroOrderThree(initialS,mu,sigma,N_approx,t0,T)
%TalayTubaroOrderThree.m: Function taking inputs defined for the problem
%along with a number of time-steps, which calls randomWalkGBM.m three times
%to obtain estimates before outputting the improved estimate with third
%order accuracy.

u_h1=randomWalkGBM(initialS,mu,sigma,N_approx,t0,T); %get first approximation by weak Euler method 
u_h2=randomWalkGBM(initialS,mu,sigma,N_approx*2,t0,T); %get second approximation with twice number of time steps as above
u_h3=randomWalkGBM(initialS,mu,sigma,N_approx*4,t0,T); %get third approximation with twice number of time steps as above
h1_approx=u_h1(length(u_h1)); %get the approximations of process at final time
h2_approx=u_h2(length(u_h2)); 
h3_approx=u_h3(length(u_h3));
extrapOrderThree=(1/3)*(8*h3_approx-6*h2_approx+h1_approx);




