clc;
close all;
clear;
%Setting up the field
nx=41;
dx=2/(nx-1);
nt=25;
dt=0.025;
c=1;
x=(0:dx:2);
%Setting up the the initial condition
u=ones(1,nx);%u_zero
%Below loop is same as: u(int64(.5/dx):int64(1/dx+1))=2;
for i=1:nx
    if (0.5<=x(i)) && (x(i)<=1)
        u(i)=2;
    end
end
%Setting up go-between parameter un
un=ones(1,nx);
%Iterating
for n=1:nt
    un=u;%setting u_zero to un
    for i=2:nx
    u(i)=un(i)-c*(dt/dx)*(un(i)-un(i-1));
    end
end
