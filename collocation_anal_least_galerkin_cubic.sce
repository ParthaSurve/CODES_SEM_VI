// 
// AEROSPACE STRUCTURES II
// Assignment 1
//--------------------------------------------------------//
// Created : 3-02-2017
//  Weighted Residual Methods : Comparing Least-Square_Collocation_Galerkin_analytical
// Ref: P.Seshu 
// Author: Partha Surve (SC14B036, Aerospace Engineering 3rd Year, IIST)


//analytical solution 
clc;
clear;
clf;


y_ana=zeros(1,21)//analytical solution
y1=zeros(1,21)
y2=zeros(1,21)
y3=zeros(1,21)
y4=zeros(1,21)

x_num=zeros(1,21)
 
//The location of the discrete points where the residual is put to zero   
    for i= 0.7//0.6:0.1:0.9
        
        x1 = 0+i
        x2 = 2-i
        
        A=[(x1*x1 -2*x1 +2),(x1*x1*x1 -2*x1*x1 +6*x1 -4);(x2*x2 -2*x2 +2),(x2*x2*x2 -2*x2*x2 +6*x2 -4)]
       
        c=[1.5*x1;1.5*x2]
       
        [x0,nsA]=linsolve(A,c)
      for x= 0:0.1:2
//Approx solution using Collocation Method        
        y2(x*10+1)= 2.5.*x +x0(1).*(x-2).*x + x0(2).*x.*x.*(x-2)       
        x=[0:0.1:2]
        
//Analytical solution of the ODE
        y1(x*10+1)=(3/sin(2)).*sin(x) + x
        
       end
        plot(x,y2,"r");
        plot(x,y1,"g");
        

end
//Least Square Method 
//weights of the Approximate solution
//W1 = 0.804165 
//W2 = 0.26727
for x= 0:0.1:2
y3(x*10+1)= 2.5.*x -0.8041.*(x-2).*x - 0.26727.*x.*x.*(x-2) 
//Galerkin method a2 = 0.97368 a3 = 0.2763 
y4(x*10+1)= 2.5.*x -0.97368.*(x-2).*x - 0.2763.*x.*x.*(x-2) 
end
x=[0:0.1:2]

plot(x,y3,"b");
plot(x,y4,"y");

title('Collocation Method', "fontsize", 5);
xlabel("Y value", "fontsize", 5);
ylabel("X value","fontsize", 5);
hl=legend(['Collocation Method';'Analytical';'LeastSquare Method']);
 
