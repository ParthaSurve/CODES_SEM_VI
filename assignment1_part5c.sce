// Assigment 1 Computational Methods for compressible flow
//--------------------------------------------------------//
// Submition date : 27-01-2017
// Contour Plot of the error between Analytical solution using series method and Explicit Forward Time Central Space scheme (FTCS) soluton for solving PDEs. Diffusion Coefficient is kept constant 
// Author: Partha Surve (SC14B036, Aerospace Engineering 3rd Year, IIST)

clc;
clear;
clf;
//initializing the parameters
L = 1  //length of the domain
t_max = 1 //end time
n_x = 8 // No of nodes
n_t = 101 // time steps
dx = L/(n_x-1); //spacial step size
dt = t_max/(n_t-1);//temporal step size

Dx=dx//fixed dx

alpha=0.001//diffusion term is varying from 0.001 to 0.1
// don't increase the value beyond 0.001 as the FTCS becomes unstable for r>0.5 

dt=0.01// keeping dt constant 


T=zeros(n_t,n_x);//initializing an array 


//Analytical Solution//
//=========================================================//
//the 'n' value can be varied to increase the accuracy 

for i = 1:n_t-1//varying time step
    for j=1:1:50//varying alpha
       T_ana(j,i)=0
    for n = 1:1:13// tolerance value of 10e-08
    U(n)= (40*(1-(-1)^n)/(n*%pi))*sin(n*%pi*Dx/L)*%e^(-alpha*dt*i*j*(n*%pi/L)^2)
    T_ana(j,i)=T_ana(j,i)+U(n)
    end    
   end
end

//===========================================================//
//Numerical Solution using Explicit Method
//============================================================//

//Boundary condition
T(1,:) = 20;
T(:,1) = 0; //dirichlet at the tip and end node
T(:,n_x) = 0;
T_ana(1,:) = 20;
T_ana(:,1) = 0; 
T_ana(:,n_x) = 0;


//initial condition at all the elements


//solving using FTCS
for j = 1:1:100
    for i = 2:1:n_x-1
        T(j+1,i) = T(j,i)+(alpha*j*dt/(dx*dx)).*(T(j,i-1)-2.*T(j,i)+T(j,i+1)); 
    
    end
end












