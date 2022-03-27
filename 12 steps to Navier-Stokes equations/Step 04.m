%Step 4 : Burgers' Equation in 1D
clc
close all
clear

%Setting up the field
nx=101;
nt=100;
nu0=0.07;
sigma=0.2;
dx=(2*pi)/(nx-1);
dt=dx*nu0;
x0=(0:dx:2*pi);
t0=0;

%Defining u_zero (sawtooth function)
syms x nu t
phi=exp(-(x-4*t)^2/(4*nu*(t+1)))+exp(-(x-4*t-2*pi)^2/(4*nu*(t+1)));
%pretty(phi);
phiprime=diff(phi,x);

u=-(2*nu*phiprime)/phi+4;
ufunc=symfun(u,[t x nu]);
u=double(ufunc(t0,x0,nu0));   %u_zero (initial time)
%plotting u_zero
%plot(x0,u,'.-','linewidth',1.1);
%axis([0 2*pi 0 10]);

%Defining go-between parameter un
un=zeros(1,nx);

%Periodic boundary condition
for n=1:nt
    un=u;
    for i=2:nx-1  %We have yet to identify u(1) and u(nx) because they are out of iteration 
        u(i)= un(i)+nu0*(dt/(dx)^2)*(un(i+1)-2*un(i)+un(i-1))...
            -(un(i))*(dt/dx)*(un(i)-un(i-1));

        u(1)=un(1)+nu0*(dt/(dx)^2)*(un(2)-2*un(1)+un(nx-1))...
            -(un(1))*(dt/dx)*(un(1)-un(nx-1)); % in this we used u(nx-1)=u(0)
        u(nx)=u(1);
    end
end

u_analytical=double(ufunc(nt*dt,x0,nu0));   %u_analytical (final time)

%Output illustration 
plot(x0,u,'.-',x0,u_analytical);
axis([0 2*pi 0 10]);
legend('Computational', 'Analytical');
