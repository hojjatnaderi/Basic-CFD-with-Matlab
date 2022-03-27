%Step 7: 2D Diffusion
clc
close all
clear
tic

%Setting up the field
nx=31;
ny=31;
nt=17;
nu=0.05;
dx=2/(nx-1);
dy=2/(ny-1);
sigma=0.25;
dt=sigma*dx*(dy/nu);

%Setting up initial conditions
u=ones(ny,nx);
u(int64(.5/dy):int64(1/dy+1),int64(.5/dx):int64(1/dx+1))=2;

%Setting up go-between parameter un
un=ones(ny,nx);

%Iteration
for n=1:nt
    un=u;
    u(2:end-1,2:end-1)=un(2:end-1,2:end-1)...
        +nu*(dt/(dx)^2)*(un(2:end-1,3:end)-2*un(2:end-1,2:end-1)+un(2:end-1,1:end-2))...
        +nu*(dt/(dy)^2)*(un(3:end,2:end-1)-2*un(2:end-1,2:end-1)+un(1:end-2,2:end-1));
    
    u(1,:)=1;
    u(end,:)=1;
    u(:,1)=1;
    u(:,end)=1;
end

%Output illustration
x=(0:dx:2);
y=(0:dy:2);
surf(x,y,u);
toc
