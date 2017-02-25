// 
// SPACE FLIGHT MECHANICS
// Assignment 1
//==================================================================//
// Created : 4-02-2017
//  
// Ref: Notes
// Author: Partha Surve (SC14B036, Aerospace Engineering 3rd Year, IIST)
// parthasurve1@gmail.com
//
clc
clear
clf
//-------------------------------------------------//
//  Constants  //
H_p = 500 //Perigee Altitude
H_a = 25000 //Apogee Altitude 
R_e = 6378 //Radius of earth
R_p = H_p + R_e // Perigee radius
R_a = H_a + R_e //Apogee Radius
a = (R_p + R_a)/2 //Semi-Major Axis
G_earth = 398600 //gravitational constant of Earth (mu)
e=(R_a-R_p)/(R_a+R_p)//eccentricity of the orbit
pi = 3.14
p=2*pi*(a^1.5)/(G_earth^0.5);//time period
n=(G_earth/(a^3))^0.5;
tolerance = 1.e-2
//--------------------------------------------------//
//Mean Anomaly Vs Time

t= 1:1:p
M=n.*t

//Evaluating Eccentric anomaly based on Mean Anomaly using Newton Rhapson
Residue = 1
for t = 1:1:p
E(t)=M(t)
while abs(Residue) > tolerance

E(t) = E(t)- (E(t) - e.*sin(E(t)) - M(t))/(1 - e.*cos(E(t)))
Residue = E(t)-e.*sin(E(t))-M(t)

end
Residue=1

//Evaluating True Anomaly from Eccentric Anomaly
v(t)=atan(((1-e*e)^.5)*sin(E(t))/(1-e*cos(E(t))),(cos(E(t))-e)/(1-e*cos(E(t))));
       if(v(t)<=0)
           v(t)=v(t)+2*pi;
       end
end
//plotting the anomalies vs time
t= 1:1:p
figure(1)
plot(t,M,'b')
plot(t,E,'r')
plot(t,v,'g')
xlabel('Time ', "fontsize", 5)
ylabel('Anomaly(Degs)', "fontsize", 5)
title('True anomaly, Eccentric anomaly, Mean anomaly over one revolution', "fontsize", 3)
hl=legend(['Mean Anomaly';'Eccentric Anomaly';'True Anomaly'],2);
//--------------------------------------------------//

//plotting radial distance

r=zeros(360,1) //radial distance


  
for V = 1:1:360 // V is true anomaly

r(V)=(a*(1-e^2)./(1+e.*cos(V*pi/180))); //Evaluating Radial Distance

end

V=  0:0.01:2*%pi;

figure(2)
polarplot(V,(a*(1-e^2)./(1+e.*cos(V)))); // Orbit of the Space craft
title('Polar plot of the Orbit ', "fontsize", 5)

//PLotting the Variation of radial distance with the true anomaly 

 V = 1:1:360 
 figure(3)
 plot(V,r,'r')
xlabel('True Anomaly(degree)', "fontsize", 5)
ylabel('Radial Distance(km)', "fontsize", 5)
title('Radial distance evolution over one revolution', "fontsize", 5)

//--------------------------------------------------------------//
// tangential Velocity
V_t = zeros(360,1) 

for V = 1:1:360 // V is true anomaly
V_t(V) = sqrt(G_earth*(2*(1+e.*cos(V*pi/180))./(a*(1-e*e)) - 1/a))

end
 V = 1:1:360 
 
 figure(4)
 plot(V,V_t,'r')
 
 xlabel('True Anomaly(degree)', "fontsize", 5)
ylabel('Tangential Velocity(km/s)', "fontsize", 5)
title('Tangential Velocity evolution over one revolution', "fontsize", 5)

//--------------------------------------------------------------//
//Flight Angle 
E = zeros(360,1)
flight_angle = zeros(360,1)
for V = 1:1:360
  A = sin(V*pi/180).*sqrt(1 - e^2)./(1 + e * cos(V*pi/180));
  B = (e + cos(V*pi/180))./(1 + e .* cos(V*pi/180));
  E(V) = atan(A/B);
  flight_angle(V) = atan((e.*sin(V*pi/180))/(1+e.*cos(V*pi/180)))
end
V = 1:1:360 
 
 figure(5)
 plot(V,flight_angle,'r')
 xlabel('True Anomaly(degree)', "fontsize", 5)
ylabel('Flight angle', "fontsize", 5)
title('', "fontsize", 5)
//==========================================================//









