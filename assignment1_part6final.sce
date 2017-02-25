// Assigment 1 Computational Methods for compressible flow
//--------------------------------------------------------//
// Submition date : 27-01-2017
// Grid convergence study
// Author: Partha Surve (SC14B036, Aerospace Engineering 3rd Year, IIST)

clc;
clear;
clf;
//initializing the parameters
L = 1  //length of the domain
t_max = 1 //end time

n_x1 = 32 // vary the No of nodes from 3 to 34 for grid convergence study

n_t = 2001 // time steps remain the same 

dx1 = L/(n_x1-1); //spacial step size

dt = t_max/(n_t-1);//temporal step size
alpha = 1; //diffusion term
t=0.01

r1=(alpha*dt/(dx1*dx1));


T1=zeros(n_t,n_x1);//initializing an array 


//Analytical Solution//
//=========================================================//
//the 'n' value can be varied to increase the accuracy 
//evaluate the analytical solution for various spatial step size
for j=2:1:n_t
    for i = 1:n_x1
   T_ana1(j,i)=0
    for n = 1:1:13// tolerance value of 10e-08
    U1(n)= (40*(1-(-1)^n)/(n*%pi))*sin(n*%pi*(i-1)*dx1/L)*%e^(-alpha*(j-1)*dt*(n*%pi/L)^2)
    T_ana1(j,i)=T_ana1(j,i)+U1(n)
    end    
   end 
end

//===========================================================//
//Numerical Solution using Explicit Method
//============================================================//



//Boundary condition
T1(1,:) = 20;
T1(:,1) = 0; //dirichlet at the tip and end node
T1(:,n_x1) = 0;
T_ana1(1,:) = 20;
T_ana1(:,1) = 0; 
T_ana1(:,n_x1) = 0;



//initial condition at all the elements


//solving using FTCS
for j = 1:1:n_t-1 
    for i = 2:1:n_x1-1
        T1(j+1,i) = T1(j,i)+r1.*(T1(j,i-1)-2.*T1(j,i)+T1(j,i+1)); 
    
    end

end

//================================================================//
// Evaluating and ploting the Error contour
//================================================================//
//constant diffusion coefficient
for j = 1:1:n_t
    
    for i = 1:1:n_x1
        err1(j,i)=T_ana1(j,i)-T1(j,i)
    end
end
  
//Calculating the RMS error
for j = 1:1:n_t
    err=0
    for i = 1:1:n_x1
        err=(T_ana1(j,i)-T1(j,i))^2 +err
    end
    RMS_err(j)=sqrt(err/n_x1)
        
end

i = [1:1:n_x1]


  
            plot(i,err1(25,:),"r");
               
title('The steady state error', "fontsize", 5);
xlabel("Nodes", "fontsize", 5);
ylabel("Error","fontsize", 5);









