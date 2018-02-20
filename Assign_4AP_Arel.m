%Assignment 4AP
%Micah Arel
%2/19/2018
%Assignment 4AP
%Solve 1D Helmholtz Equation using Tri-Diagonal Method
%2nd Derivative approximation based on Center Difference Formula
%Total of N nodes between 0 and L; x=0 N=0, x=L N+1
clear
clc
L = 1; %Right endpoint
V0=1;
v=1;
A=1; %Given
k = 10; %Coefficient of u(x) in ODE
BC01 = V0; %Boundary condition node 0 (left side)for problem 1
BC02 = v; %Boundary condition node 0 (left side) for problem 2
BCL = 0; %Boundary condition node N+1 (right side) for both problems
N=89;
h = L/(N+1); %Step size
D = -(2+k*k*h*h); %Coefficient

%Create 4 vectors of length N, for diagonals and RHS of system
MD1 = D.*ones(1,N); %Main Diagonal
LD1 = ones(1,N); %Lower Diagonal
UD1 = ones(1,N); %Upper Diagonal
RHS1 = (h^2)*ones(N,1); %RHS vector (given f(x))
RHS1(1) = RHS1(1) - BC01; %Input boundary conditions to Right Hand Side
RHS1(N) = RHS1(N) - BCL;
%Create Problem 1 upper triangle forward elimination
for i = 2:(N)
 MD1(i) = MD1(i) - ( (LD1(i)/MD1(i-1)) * UD1(i-1) );
 RHS1(i) = RHS1(i) - ( (LD1(i)/MD1(i-1)) * RHS1(i-1) );
end
%Solve Problem 1 system backwards elimination
U1(N)=RHS1(N)/MD1(N);
for j = N-1:-1:1
 U1(j) = ( RHS1(j) - UD1(j)*U1(j+1) ) / MD1(j);
end
%Complete Problem 1 U vector
U1 = [BC01 U1 BCL];

%Display to if N and 2N are approximately =
Z1 = U1.';
X1=zeros(N,1);
for j = 2:(N+2)
    X1(j)=X1(j-1)+h;
end
Y1=[X1 Z1]; %Checks accuracy at N and at 2N

%Loop used to evaluate the difference of the current and previous
%solution at the previous solutions spatial location
%Analytical solution
x1=0:h:L;
Ua=( (1 - ( ( sinh(k*(L-x1)) + sinh(k*x1) ) / sinh(k*L) )) * (A/k^2) ) + BC01*(sinh(k*(L-x1))/sinh(k*L));

%Plot of approximation and exact solutions
x=0:h:L;
figure(1)
plot(x1,Ua,x,U1,'*')
lgd=legend('Exact Solution','Approximate Solution');
title(lgd,{'\Delta = L / (N+1) = 0.0111';'N = 89';'k = 10'})
title({'1-D Helmholtz Equation';'Dirichlet BC (Problem 1)'})
xlabel('X')
ylabel('U')

%Problem 2

%Create 4 vectors of length N, for diagonals and RHS of system
MD2 = D.*ones(1,N); %Main Diagonal
LD2 = ones(1,N); %Lower Diagonal
UD2 = ones(1,N); %Upper Diagonal
RHS2 = (h^2)*ones(N,1); %RHS vector (given f(x))
RHS2(1) = RHS2(1) + 2*h*BC02; %Input boundary conditions to RHS
RHS2(N) = RHS2(N) - BCL;

%Create Problem 2 upper triangle forward elimination
for i = 2:(N)
 MD2(i) = MD2(i) - ( (LD2(i)/MD2(i-1)) * UD2(i-1) );
 RHS2(i) = RHS2(i) - ( (LD2(i)/MD2(i-1)) * RHS2(i-1) );
end
%Solve Problem 2 system backwards elimination
U2(N)=RHS2(N)/MD1(N);
for j = N-1:-1:1
 U2(j) = ( RHS2(j) - UD2(j)*U2(j+1) ) / MD2(j);
end
%Complete Problem 1 U vector
U2 = [BC02 U2 BCL];

%Loop used to evaluate the difference of the current and previous
%solution at the previous solutions spatial location
%Analytical solution
x2=0:h:L;
Ub=( (1 - (  cosh(k*(L-x2)) / cosh(k*L) )) * (A/k^2) ) + BC02/k*(sinh(k*(L-x2))/cosh(k*L));

%Plot of approximation and exact solutions
x=0:h:L;
figure(2)
plot(x1,Ub,x,U2,'*')
lgd=legend('Exact Solution','Approximate Solution');
title(lgd,{'\Delta = L / (N+1) = 0.0111';'N = 89';'k = 10'})
title({'1-D Helmholtz Equation';'Neumann BC (Problem 2)'})
xlabel('X')
ylabel('U')
