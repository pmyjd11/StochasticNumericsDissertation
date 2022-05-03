%Script which produces a comparison of the scaled random walk with Brownian
%motion (figure 16 in the report, page 98).

%Form 3 time vectors of differing size
t1=linspace(0,1,100);
h1=t1(2)-t1(1); 
t2=linspace(0,1,500);
h2=t2(2)-t2(1); 
t3=linspace(0,1,1000);
h3=t3(2)-t3(1); 
%Form Brownian motion on finest mesh
B=zeros(1,length(t3));
B(1)=0; %Brownian motion starts from zero
for i=1:length(t3)-1
    B(i+1)=B(i)+normrnd(0,sqrt(h3)); %each increment is N(0,h)
end
%Form scaled random walk
nu1=zeros(1,length(t1));
for i=1:length(t1)-1
    J=rand;
    if J<0.5
        nu1(i+1)=nu1(i)-1; %subtract 1
    else
        nu1(i+1)=nu1(i)+1; %add 1
    end
end
nu1=nu1.*sqrt(h1); %scale the walk
%Repeat for the other two
nu2=zeros(1,length(t2));
for i=1:length(t2)-1
    J=rand;
    if J<0.5
        nu2(i+1)=nu2(i)-1; %subtract 1
    else
        nu2(i+1)=nu2(i)+1; %add 1
    end
end
nu2=nu2.*sqrt(h2);
nu3=zeros(1,length(t3));
for i=1:length(t3)-1
    J=rand;
    if J<0.5
        nu3(i+1)=nu3(i)-1; %subtract 1
    else
        nu3(i+1)=nu3(i)+1; %add 1
    end
end
nu3=nu3.*sqrt(h3);
%Now plot
plot(t3,B,'k');
hold on
plot(t1,nu1,'g');
plot(t2,nu2,'r');
plot(t3,nu3,'b');
xlabel('t');
ylabel('Process at time t');
title('Comparison of scaled random walk with Brownian motion');
legend('Brownian motion (1000 points)','Scaled walk (100 points)','Scaled walk (500 points)','Scaled walk (1000 points)');
