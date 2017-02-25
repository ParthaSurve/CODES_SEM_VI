// 
// AEROSPACE STRUCTURES II
// Assignment 1
//--------------------------------------------------------//
// Created : 3-02-2017
//  Weighted Residual Methods : Comparing Least-Square_Collocation_Galerkin_analytical
// Ref: P.Seshu 
// Author: Partha Surve (SC14B036, Aerospace Engineering 3rd Year, IIST)

//====================================================//
//Constants//
clear;
clc;

L = 50 ; //length of the beam 
E = 100; //youngs miodulus
P = 100; //axial load at the tip

n= 4 // nodes
elem =n-1;//elements
dx= L/elem;//length of each element 

//distributed force
function y=q(x)
    
    y=x^2
   
endfunction

// LOCAL STIFFNESS MATRIX

K = zeros(elem,2,2)

for i = 1:1:elem
    
    x= (i-1)*dx
    r = -0.1 * x + 10; // equation of the slope of the tapered bar
    A = %pi * r* r     // Area of the cross section 
        
  K(i,2,2) = A*E/dx
  K(i,1,1) = A*E/dx
  K(i,2,1) = -A*E/dx
  K(i,1,2) = -A*E/dx
    
end


//GLOBAL STIFFNESS MATRIX

K_stiff = zeros(n,n)

    for i= 1:1:elem

        
                  if (i< elem)
                      a = K(i,1,1)
                      b = K(i,1,2)
                     
                     
                     K_stiff(i,i) = a + K_stiff(i,i) 
                     K_stiff(1+i,1+i) = a + K_stiff(1+i,1+i)
                     
                     K_stiff(i,i+1) = b  
                     K_stiff(i+1,i) = b  
                     
                     
                  end
                 
                  if (i==elem)
                      
                      c = K(i,1,2)
                      d = K(elem,2,2)
                     
                     K_stiff(i,i) = d + K_stiff(i,i) 
                     K_stiff(1+i,1+i) = d + K_stiff(1+i,1+i)
                     
                     K_stiff(i,i+1) = c  
                     K_stiff(i+1,i) = c 
                  end
                     
 end




// applying the boundary condition at the the first node . u1 = 0 .
// The first column and the first row of the global stiffness gets removed 
// hence we get a reduced stiffness matrix

//REDUCED FORCE MATRIX

F =zeros(n-1,1)
F(n-1,1) = P

//DISTRIBUTED FORCE

N = zeros(elem,2)

for i =1:1:n-1
        
    x= (i-1)*dx
    r = -0.1 * x + 10; // equation of the slope of the tapered bar
    A = %pi * r* r     // Area of the cross section 
     x1=0 +  (i-1)*dx
     x2=  i*dx  
      //2 POINT GAUSS-QUADRATURE 
      N(i,1) = (A*dx/2)*(0.78665*q(0.5*(1.57733*x1 + 0.4227* x2))+(0.21135*q(0.5*(0.4227*x1 + 1.57733* x2))))
      
      N(i,2) = (A*dx/2)*(0.21135*q(0.5*(1.57733*x1 + 0.4227* x2))+(0.78665*q(0.5*(0.4227*x1 + 1.57733* x2))))
      
      //refer to //http://www.engr.uvic.ca/~mech420/Lecture20.pdf
    //for the gauss qudrature formula used
    end

Q= zeros(n,1)
for j =2:1:n-1
        
         Q(1) =  N(1,1)
         Q(n) =  N(elem,2)
      Q(j) = N(j-1,2) + N(j,1) 
    
    
    end


//TOTAL FORCE ACTING ON THE ROD
F_total =zeros(n-1,1)
for j =1:1:n-1
        
        F_total(j) = Q(j+1) 
        F_total(n-1) = F_total(n-1) + P
    
    end



//REDUCED STIFFNES MATRIX
 K_reduced=zeros(n-1,n-1)
for i=1:1:n-1
    for j =1:1:n-1
        
         
      K_reduced(i,j) = K_stiff(i+1,j+1)
    
    
    end
        
end


// for evaluating th displacements at each of the nodes

[u,nsA]=linsolve(K_reduced,-F_total)




//REACTION FORCE


R =0 //reaction force at the fixed boundary condtion


for j =1:1:n-1

    R = R + u(j)* K_stiff(1,j+1) 
    
end

R = R - Q(1,1)








