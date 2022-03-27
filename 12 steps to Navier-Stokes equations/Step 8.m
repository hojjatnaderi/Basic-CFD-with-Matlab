%Step 8: Burgers' Equation in 2D
clc
close all
clear
tic
%
nx=41;
ny=41;
nt=120;
dx=2/(nx-1);
dy=2/(ny-1);
sigma=0.0009;
nu=0.01;
dt=sigma*dx*(dy/nu);
%
u=ones(ny,nx);
v=ones(ny,nx);
u(int64(.5/dy):int64(1/dy+1),int64(.5/dx):int64(1/dx+1))=2;
v(int64(.5/dy):int64(1/dy+1),int64(.5/dx):int64(1/dx+1))=2;
%
for n=1:nt
    un=u;
    vn=v;
    for j=2:ny-1
        for i=2:nx-1
            u(j,i)=un(j,i)-(dt/dx)*un(j,i)*(un(j,i)-un(j,i-1))...
                -(dt/dy)*vn(j,i)*(un(j,i)-un(j-1,i))...
                +nu*(dt/(dx)^2)*(un(j,i+1)-2*un(j,i)+un(j,i-1))...
                +nu*(dt/(dy)^2)*(un(j+1,i)-2*un(j,i)+un(j-1,i));
            
            v(j,i)=vn(j,i)-(dt/dx)*un(j,i)*(vn(j,i)-vn(j,i-1))...
                -(dt/dy)*vn(j,i)*(vn(j,i)-vn(j-1,i))...
                +nu*(dt/(dx)^2)*(vn(j,i+1)-2*vn(j,i)+vn(j,i-1))...
                +nu*(dt/(dy)^2)*(vn(j+1,i)-2*vn(j,i)+vn(j-1,i));
            
            u(1,:)=1;
            u(end,:)=1;
            u(:,1)=1;
            u(:,end)=1;
            
            v(1,:)=1;
            v(end,:)=1;
            v(:,1)=1;
            v(:,end)=1;
        end
    end
end
%
x=(0:dx:2);
y=(0:dy:2);
surf(x,y,u);
toc
