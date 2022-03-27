clc
close all
clear
%Setting up the field
nx=41;
nt=20;
dx=(2-0)/(nx-1);
dt=0.025;
%Setting up u_zero
u=ones(1,nx);
u(int64((0.5-0)/dx): int64((1-0)/dx+1))=2;
%Plotting u_zero
u_zero=u;
x=(0:dx:2);
%plot(x,u_zero)
%hold on
%Defining un parameter
un=ones(1,nx);
%Iterating
 for n=1:nt
     un=u;
     for i=2:nx
         u(i)=un(i)-(un(i))*(dt/dx)*(un(i)-un(i-1));
     end
 end
 %Plotting output
 plot(x,u)
