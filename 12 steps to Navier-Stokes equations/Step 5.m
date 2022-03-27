%Step 5: 2-D Linear Convection (using "array operations" instead of "loops")
clc
close all
clear
tic %this is for runtime measuring

%Setting up the field
nx=81;
ny=81;
nt=100;
c=1;
dx=(2-0)/(nx-1);
dy=(2-0)/(ny-1);
sigma=0.2;
dt=sigma*dx;

%Initial conditions
u=ones(ny,nx);
u(int64(0.5/dy):int64(1/dy+1),int64(0.5/dx):int64(1/dx+1))=2;
uzero=u;

%Plotting initial conditions
x=(0:dx:2);
y=(0:dy:2);
surf(x,y,uzero);

%un is go-between parameter
un=ones(ny,nx); 

%2nd method: Iterating in two dimensions by array operations (1st method is done by for-loops)
for n=1:nt
    un=u;
    u(2:end,2:end)=un(2:end,2:end)-(c*dt/dx)*(un(2:end,2:end)-un(2:end,1:end-1))...
        -(c*dt/dy)*(un(2:end,2:end)-un(1:end-1,2:end));

    u(1,:)=1;
    u(end,:)=1;
    u(:,1)=1;
    u(:,end)=1;
     
end

%Output
surf(x,y,u);
toc %this is for runtime measuring (you can notice considerable reduction in runtime when not using loops)
