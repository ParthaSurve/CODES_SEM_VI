clear;
clc
clf

pi = 2*3.14
M=2.9288 //mean anomaly angle
e = 0.9 //eccentricity
tolerance = 1.e-6 //error tolerance

//Select a starting value for E:


//Iterate using the Newton Rhapson Method

Residue = 1
E= ones(360,1)
E1= ones(360,1)
for M = 1:1:360

while abs(Residue) > tolerance

E(M) = E(M)- (E(M) - e.*sin(E(M)) - (M*pi/360))/(1 - e.*cos(E(M)))
Residue = E(M)-e.*sin(E(M))-M*pi/360

end
//----------------------//
//Series Solution 
//Ref : Curtis
E1(M)= M*pi/360 + e*sin(M*pi/360) + ((e^2)/2)*sin(2*M*pi/360) + ((e^3)/8)*(3*sin(3*M*pi/360)-sin(M*pi/360))
//------------------------//
E(M)=E(M)*360/pi;
E1(M)=E1(M)*360/pi;
Residue =1

end
M=[1:1:360]
        plot(M,E,"r");
        plot(M,E1,"b");
        title('Solution To Kepler Equation', "fontsize", 5);
xlabel("Eccentric Anomaly (degs)", "fontsize", 5);
ylabel("Mean Anomaly (degs)","fontsize", 5);
hl=legend(['Numerical Solution';'Series Solutoin'],2);
