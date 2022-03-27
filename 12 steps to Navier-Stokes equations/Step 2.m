clc
clear all
close all
nx=41;
nt=20;
dx=(2-0)/(nx-1);
dt=0.025;
%defining u_zero
u=ones(1,nx);
u(int64((0.5-0)/dx): int64((1-0)/dx+1))=2;

u_zero=u;
x=(0:dx:2);
%plot(x,u_zero)
%hold on
%defining un middle parameter
un=ones(1,nx);
 for n=1:nt
     un=u;
     for i=2:nx
         u(i)=un(i)-(un(i))*(dt/dx)*(un(i)-un(i-1));
     end
 end
 
 plot(x,u)
