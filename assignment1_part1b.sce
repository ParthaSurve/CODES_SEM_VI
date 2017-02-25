// Assigment 1 Computational Methods for compressible flow
//--------------------------------------------------------//
// January 17-01-2017
// Explicit Forward Time Central Space scheme (FTSC) for solving PDEs
// Author: Partha Surve (SC14B036, Aerospace Engineering, IIST)


clc;
clear;
clf;
//==========================//
//initializing the parameters
//==========================//

L = 1  //length of the domain
t_max = 1 //end time

n_x = 8 // No of nodes
n_t = 101 // time steps

dx = L/(n_x-1); //spacial step size
dt = t_max/(n_t-1);//temporal step size

alpha = 1; //diffusion term
r=(alpha*dt/(dx*dx));
T=zeros(n_t,n_x);//initializing an array 

//Boundary condition
T(:,1) = 0; //dirichlet at the tip and end node




for k = 2:1:n_x
    T(1,k)=6*sin(%pi*(k-1)*dx/L);
end

//solving using FTCS
for j = 1:1:n_t-1 
    for i = 2:1:n_x-1
        T(j+1,i) = T(j,i)+r.*(T(j,i-1)-2.*T(j,i)+T(j,i+1));
        T(j+1,n_x)= T(j+1,n_x-1)-dx 
    
    end
end
    
//Plotting for each time step    
i=[1:1:n_x]
for j = 1:1:n_t   
       plot(i,T(j,:),"b");
 end    
title('1-D Heat Eqn-Explicit-Neumann condition', "fontsize", 5);
xlabel("Nodes", "fontsize", 5);
ylabel("funciton value","fontsize", 5);





