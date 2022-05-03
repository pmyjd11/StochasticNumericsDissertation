%Code that simulates one-dimensional Brownian motion, providing three
%magnifications of a trajectory over the domain [0,1].
rng(1);
t=linspace(0,1,100000); %discretise the domain [0,1]
tZoom=linspace(0.5,0.6,10000);
tZoomPlus=linspace(0.55,0.552,200);
h=t(2)-t(1); %the distance between any two points
hZoom=tZoom(2)-tZoom(1);
hZoomPlus=tZoomPlus(2)-tZoomPlus(1);
B=zeros(1,length(t)); 
BZoom=zeros(1,length(tZoom));
BZoomPlus=zeros(1,length(tZoomPlus));
B(1)=0; %Brownian motion starts from zero
for i=1:length(t)-1
    B(i+1)=B(i)+normrnd(0,sqrt(h)); %each increment is N(0,h)
end
for j=1:10000
    BZoom(j)=B((length(t)/2)+j-1);
end
for k=1:200
    BZoomPlus(k)=B(55000+k-1);
end
subplot(3,1,1);
plot(t,B);
xlabel('t');
ylabel('B(t)');
title('Simulation of Brownian motion on [0,1]');
subplot(3,1,2);
plot(tZoom,BZoom);
xlim([0.5,0.6]);
xlabel('t');
ylabel('B(t)');
title('Simulation of Brownian motion on [0.5,0.6]');
subplot(3,1,3);
plot(tZoomPlus,BZoomPlus);
xlim([0.55,0.552]);
xlabel('t');
ylabel('B(t)');
title('Simulation of Brownian motion on [0.55,0.552]');
