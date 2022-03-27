%Step 3: Diffusion Equation in 1-D
clc
close all
clear
%Setting up the field
nx=41;
nt=20;
nu=0.3;
sigma=0.2;
dx=(2-0)/(nx-1);
dt=sigma*(dx^2/nu);
%Defining hat function for u_zero
u=ones(1,nx);
u(int64((0.5-0)/dx):int64((1-0)/dx+1))=2;
%Plotting u_zero
u_zero=u;
x=(0:dx:2);
plot(x,u_zero)
hold on
%Defining go-between parameter un
un=ones(1,nx);
%Iterating through time and space

for n=1:nt
    un=u;
    for i=2:nx-1
        u(i)= un(i)+nu*(dt/(dx)^2)*(un(i+1)-2*un(i)+un(i-1));
    end
end
%Output
plot(x,u)
