function randomWalkApprox = randomWalkGBM(initialS,mu,sigma,N_approx,t0,T)
%randomWalkGBM: Function that takes inputs from the problem plus a number
%of time-steps, and outputs the weak Euler approximation. 

hApprox=(T-t0)/N_approx; %find the size of each time-step for approximate 
randomWalkApprox=zeros(1,N_approx+1); %initialise the output vectors
randomWalkApprox(1)=initialS; %assign the first value in each vector to initialS
for i=1:N_approx
    if rand>0.5
    randomWalkApprox(i+1)=randomWalkApprox(i)+hApprox*mu*randomWalkApprox(i)+sigma*randomWalkApprox(i)*sqrt(hApprox); %calculate next element in approximation
    else
    randomWalkApprox(i+1)=randomWalkApprox(i)+hApprox*mu*randomWalkApprox(i)-sigma*randomWalkApprox(i)*sqrt(hApprox); 
    end
end

