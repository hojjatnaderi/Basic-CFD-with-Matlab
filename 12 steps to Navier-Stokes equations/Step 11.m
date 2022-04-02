%Navier Stokes - Cavity flow
clc
close all
clear
tic
%Field parameters
nx=41;
ny=41;
nt=100;
nit=50;
dx=2/(nx-1);
dy=2/(ny-1);

x=(0:dx:2);
y=(0:dy:2);

rho=1;
nu=0.1;
dt=0.001;

%Initial conditions
u=zeros(ny,nx);
v=zeros(ny,nx);
p=zeros(ny,nx);
b=zeros(ny,nx);

un=zeros(ny,nx);
vn=zeros(ny,nx);
pn=zeros(ny,nx);

%Iteration
for n=1:nt
    
    un=u;
    vn=v;
    pn=p;
    
    u(2:end-1,2:end-1)=un(2:end-1,2:end-1)...
        -un(2:end-1,2:end-1)*(dt/dx).*(un(2:end-1,2:end-1)-un(2:end-1,1:end-2))... % .* is important but here it shows no effect than *
        -vn(2:end-1,2:end-1)*(dt/dy).*(un(2:end-1,2:end-1)-un(1:end-2,2:end-1))...
        -(dt/(rho*2*dx))*(pn(2:end-1,3:end)-pn(2:end-1,1:end-2))...
        +nu*((dt/dx^2)*(un(2:end-1,3:end)-2*un(2:end-1,2:end-1)+un(2:end-1,1:end-2))...
        +(dt/dy^2)*(un(3:end,2:end-1)-2*un(2:end-1,2:end-1)+un(1:end-2,2:end-1)));
    
    v(2:end-1,2:end-1)=vn(2:end-1,2:end-1)...
        -un(2:end-1,2:end-1)*(dt/dx).*(vn(2:end-1,2:end-1)-vn(2:end-1,1:end-2))...
        -vn(2:end-1,2:end-1)*(dt/dy).*(vn(2:end-1,2:end-1)-vn(1:end-2,2:end-1))...
        -(dt/(rho*2*dy))*(pn(3:end,2:end-1)-pn(1:end-2,2:end-1))...
        +nu*((dt/dx^2)*(vn(2:end-1,3:end)-2*vn(2:end-1,2:end-1)+vn(2:end-1,1:end-2))...
        +(dt/dy^2)*(vn(3:end,2:end-1)-2*vn(2:end-1,2:end-1)+vn(1:end-2,2:end-1)));
    
    b(2:end-1,2:end-1)=(1/(2*dt))*(((un(2:end-1,3:end)-un(2:end-1,1:end-2))/dx)+((vn(3:end,2:end-1)-vn(1:end-2,2:end-1))/dy))...
        -((un(2:end-1,3:end)-un(2:end-1,1:end-2))/(2*dx)).*((un(2:end-1,3:end)-un(2:end-1,1:end-2))/(2*dx))...
        -2*((un(3:end,2:end-1)-un(1:end-2,2:end-1))/(2*dy)).*((vn(2:end-1,3:end)-vn(2:end-1,1:end-2))/(2*dx))...
        -((vn(3:end,2:end-1)-vn(1:end-2,2:end-1))/(2*dy)).*((vn(3:end,2:end-1)-vn(1:end-2,2:end-1))/(2*dy));
    
    p(2:end-1,2:end-1)=((pn(2:end-1,3:end)+pn(2:end-1,1:end-2))*(dy^2)+(pn(3:end,2:end-1)+pn(1:end-2,2:end-1))*(dx^2))...
        /(2*(dx^2+dy^2))-(rho*(dx^2)*(dy^2)*b(2:end-1,2:end-1))/(2*(dx^2+dy^2));
    
    %Boundary conditions
    u(end,:)=1; %lid movement
    u(1,:)=0;
    v(end,:)=0;
    v(1,:)=0;
    
    u(:,1)=0;
    u(:,end)=0;
    v(:,1)=0;
    v(:,end)=0;
    
    p(end,:)=0;
    p(1,:)=p(2,:); %dp/dy
    
    p(:,1)=p(:,2); %dp/dx
    p(:,end)=p(:,end-1); %dp/dx
                                  
end
%Output
contourf(x,y,p);
hold on
quiver(x,y,u,v,3);
%streamline(x,y,u,v);
hold off
toc                 
