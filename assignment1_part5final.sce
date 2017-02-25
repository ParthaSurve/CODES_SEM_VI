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

n_x = 8 // No of nodes
n_t = 101 // time steps

dx = L/(n_x-1); //spacial step size
dt = t_max/(n_t-1);//temporal step size

alpha = 0.001; //diffusion term

t=0.01


T=zeros(100,n_t,n_x);//initializing an array 


//Analytical Solution//
//=========================================================//
//the 'n' value can be varied to increase the accuracy 

for k=1:1:100
    for j=2:1:n_t
        for i = 1:n_x
        T_ana(k,j,i)=0
            for n = 1:1:13// tolerance value of 10e-08
            U(n)= (40*(1-(-1)^n)/(n*%pi))*sin(n*%pi*(i-1)*dx/L)*%e^(-alpha*k*(            j-1)*dt*(n*%pi/L)^2)
             T_ana(k,j,i)=T_ana(k,j,i)+U(n)
            end    
        end
    end
end
//===========================================================//
//Numerical Solution using Explicit Method
//============================================================//



//Boundary condition
T(:,1,:) = 20;
T(:,:,1) = 0; //dirichlet at the tip and end node
T(:,:,n_x) = 0;
T_ana(:,1,:) = 20;
T_ana(:,:,1) = 0; 
T_ana(:,:,n_x) = 0;


//initial condition at all the elements


//solving using FTCS

for k = 1:1:100
    for j = 1:1:n_t-1 
        for i = 2:1:n_x-1
            T(k,j+1,i) = T(k,j,i)+(alpha*k*dt/(dx*dx)).*(T(k,j,i-1)-2.*T(k,j,i)+T(k,j,i+1)); 
        end
    end
end
  
//================================================================//
// Evaluating and ploting the Error contour
//================================================================//
//constant diffusion coefficient
    
for  k = 1:1:100  
    for j = 1:1:n_t
        for i = 1:1:n_x
            err(k,j,i)=T_ana(k,j,i)-T(k,j,i) 
        end
    end
end  


//keeping the diffusion term constant and varying the time and spatial step
 for j = 1:1:n_t
        for i = 1:1:n_x
            err1(j,i)=err(100,j,i)
        end
    end
//varying the diffusion term and keeping the time step at 0.5
 for k = 1:1:100
        for i = 1:1:n_x
            err2(k,i)=err(k,50,i)
        end
    end
//keeping the spatial step constant
 for k = 1:1:100
        for j = 1:1:n_t
            err3(k,j)=err(k,j,n_x/2)
        end
    end


j = [1:1:n_t]
i = [1:1:n_x]
k = [1:1:100]

// turn on the required plot by turning off the comment
//plot3d1(j,i,err1(j,i))
//plot3d1(k,i,err2(k,i))
plot3d1(k,j,err3(k,j))






