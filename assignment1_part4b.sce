// Assigment 1 Computational Methods for compressible flow
//--------------------------------------------------------//
// January 27-01-2017
// Steady state time for Implicit Backward Time Central Space scheme (BTSC) for solving PDEs
// Author: Partha Surve (SC14B036, Aerospace Engineering 3rd year, IIST)

clc;
clear;
clf;
//initializing the parameters
//--------------------------//
L = 1  //length of the domain
t_max = 1 //end time
n_x = 8 // No of nodes
n_t = 101 // time steps
dx = L/(n_x-1); //spacial step size
dt = t_max/(n_t-1);//temporal step size
alpha=1;

//initializing an array 
A=zeros(n_x,n_x);
T=zeros(n_t,n_x);
T_old=ones(1,n_x);

//components of the tridiagonal matrix
a= 1+(2*alpha*dt/(dx*dx));
b= -alpha*dt/(dx*dx);
c= -alpha*dt/(dx*dx);



//Boundary condition
//-------------------//


for i = 1:1:n_x
    A(i,i)=a;
    
end
for i = 1:1:n_x-1
    A(i,i+1)=b;
    A(i+1,i)=c;

end

//Boundary condition
A(1,1) = 1; //dirichlet at the tip and end node
A(2,1) = 0;
A(n_x,n_x) = 1;
A(n_x-1,n_x) = 0;

//temperature at the nodes
T_old(1) = 0;
T_old(n_x) = 0;

for i = 2:1:n_x-1// initial temp of rod
    T_old(i)=20;
end

//inverting the matrix to get the temp at the nodes
//TDMA can be used for more effiecient calculation.... need to work on that
for i = 1:1:n_t-1
    T(1,:)=T_old
    T(i+1,:)=T(i,:)*inv(A);
end


//Calculating the RMS error
for j = 1:1:n_t-1
    err=0
    for i = 1:1:n_x
        err=(T(j,i)-T(j+1,i))^2 +err
    end
    RMS_err(j)=sqrt(err/n_x)
    TIME=j*dt
        if RMS_err(j) < n_x*10^(-3)*dt
          break
          
          end
      
end



