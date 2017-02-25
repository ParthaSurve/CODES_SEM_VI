// Assigment 1 Computational Methods for compressible flow
//--------------------------------------------------------//
// Submition date : 27-01-2017
// Comparison between Analytical solution using series method and Explicit Forward Time Central Space scheme (FTSC) solutonfor solving PDEs
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
alpha = 1; //diffusion term

t=0.01//time at which we are evaluating the analytical funciton

r=(alpha*dt/(dx*dx));
T=zeros(n_t,n_x);//initializing an array 


//Analytical Solution//
//=========================================================//
//the 'n' value can be varied to increase the accuracy 
for i = 1:n_x
   T_ana(i)=0
    for n = 1:1:13// tolerance value of 10e-08
    U(n)= (40*(1-(-1)^n)/(n*%pi))*sin(n*%pi*(i-1)*dx/L)*%e^(-alpha*t*(n*%pi/L)^2)
    T_ana(i)=T_ana(i)+U(n)
    end    
end
//Plotting the funciton value at a particular time step for all the nodes
 
j=[1:1:n_x]
 
       plot(j,T_ana(j),"r");
     

//===========================================================//
//Numerical Solution using Explicit Method
//============================================================//



//Boundary condition
T(:,1) = 0; //dirichlet at the tip and end node
T(:,n_x) = 0;


//initial condition at all the elements
for k = 2:1:n_x-1
    T(:,k)=20;
end

//solving using FTCS
for j = 1:1:n_t-1 
    for i = 2:1:n_x-1
        T(j+1,i) = T(j,i)+r.*(T(j,i-1)-2.*T(j,i)+T(j,i+1)); 
    
    end
end
    
//Plotting for each time step    
i=[1:1:n_x]
  
       plot(i,T(2,:),"b");//plotting for t=0.01
    
       

   
title('Comparison-Analytical and Numerical solution ', "fontsize", 5);
xlabel("Nodes", "fontsize", 5);
ylabel("function value","fontsize", 5);
   
//================================================================//













