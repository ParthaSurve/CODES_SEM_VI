// Assigment 1 Computational Methods for compressible flow
//--------------------------------------------------------//
// Submition date : 27-01-2017
// Contour Plot of the error between Analytical solution using series method and Explicit Forward Time Central Space scheme (FTSC) soluton for solving PDEs. Diffusion Coefficient is kept constant 
// Author: Partha Surve (SC14B036, Aerospace Engineering 3rd Year, IIST)

clc;
clear;
clf;
//initializing the parameters
L = 1  //length of the domain
t_max = 1 //end time

n_x1 = 4 // No of nodes
n_x2 = 8
n_x3 = 16

n_t = 801 // time steps remain the same 

dx1 = L/(n_x1-1); //spacial step size
dx2 = L/(n_x2-1)
dx3 = L/(n_x3-1)

dt = t_max/(n_t-1);//temporal step size
alpha = 1; //diffusion term
t=0.01

r1=(alpha*dt/(dx1*dx1));
r2=(alpha*dt/(dx2*dx2));
r3=(alpha*dt/(dx3*dx3));

T1=zeros(n_t,n_x1);//initializing an array 
T2=zeros(n_t,n_x3)
T3=zeros(n_t,n_x3)

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
   
    for i = 1:n_x2
   T_ana2(j,i)=0
    for n = 1:1:13// tolerance value of 10e-08
    U2(n)= (40*(1-(-1)^n)/(n*%pi))*sin(n*%pi*(i-1)*dx2/L)*%e^(-alpha*(j-1)*dt*(n*%pi/L)^2)
    T_ana2(j,i)=T_ana2(j,i)+U2(n)
    end    
   end
   
    for i = 1:n_x3
   T_ana3(j,i)=0
    for n = 1:1:13// tolerance value of 10e-08
    U3(n)= (40*(1-(-1)^n)/(n*%pi))*sin(n*%pi*(i-1)*dx3/L)*%e^(-alpha*(j-1)*dt*(n*%pi/L)^2)
    T_ana3(j,i)=T_ana3(j,i)+U3(n)
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


T2(1,:) = 20;
T2(:,1) = 0; //dirichlet at the tip and end node
T2(:,n_x2) = 0;
T_ana2(1,:) = 20;
T_ana2(:,1) = 0; 
T_ana2(:,n_x2) = 0;

T3(1,:) = 20;
T3(:,1) = 0; //dirichlet at the tip and end node
T3(:,n_x3) = 0;
T_ana3(1,:) = 20;
T_ana3(:,1) = 0; 
T_ana3(:,n_x3) = 0;


//initial condition at all the elements


//solving using FTCS
for j = 1:1:n_t-1 
    for i = 2:1:n_x1-1
        T1(j+1,i) = T1(j,i)+r1.*(T1(j,i-1)-2.*T1(j,i)+T1(j,i+1)); 
    
    end

for i = 2:1:n_x2-1
        T2(j+1,i) = T2(j,i)+r2.*(T2(j,i-1)-2.*T2(j,i)+T2(j,i+1)); 
    
    end

for i = 2:1:n_x3-1
        T3(j+1,i) = T3(j,i)+r3.*(T3(j,i-1)-2.*T3(j,i)+T3(j,i+1)); 
    
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

    for i = 1:1:n_x2
        err2(j,i)=T_ana2(j,i)-T2(j,i)
    end
    
    for i = 1:1:n_x3
        err3(j,i)=T_ana3(j,i)-T3(j,i)
    end

end
  


i = [1:1:n_x1]
j = [1:1:n_x2]
k = [1:1:n_x3]


  
            plot(i,err1(25,:),"r");
            plot(j,err2(25,:),"r");
            plot(k,err3(25,:),"r");    
title('The steady state error', "fontsize", 5);
xlabel("Nodes", "fontsize", 5);
ylabel("Error","fontsize", 5);









