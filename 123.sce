syms a e; a is semimajor and e is eccentricity
[a,e] = solve(a*(1-e)==6375+300,a*(1+e)==6375+10000,a,e); % earth radius = 6375km,
hp=300, ha=10000
b=a*sqrt(1-e.^2); %b is semiminor
b=double(b);
a=double(a);
e=double(e);
u=0:pi/100:2*pi; %u is true anomaly
r=(a*(1-e^2)./(1+e.*cos(u))); %r is radial distance from earth's center
figure(1);
plot(u,r); %Radial distance evolution over one revolution
xlabel('true anomaly(degree)')
ylabel('radial distance(km)')
title('Radial distance evolution over one revolution')
