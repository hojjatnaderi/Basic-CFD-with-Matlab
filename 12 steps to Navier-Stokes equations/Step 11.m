%Navier Stokes - Cavity flow
%For discussion on how to improve this code refer to Step 06
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
    for j=2:ny-1
        for i=2:nx-1
            
            u(j,i)=un(j,i)-un(j,i)*(dt/dx)*(un(j,i)-un(j,i-1))...
                -vn(j,i)*(dt/dy)*(un(j,i)-un(j-1,i))...
                -(dt/(rho*2*dx))*(pn(j,i+1)-pn(j,i-1))...
                +nu*((dt/dx^2)*(un(j,i+1)-2*un(j,i)+un(j,i-1))...
                +(dt/dy^2)*(un(j+1,i)-2*un(j,i)+un(j-1,i)));
            
            v(j,i)=vn(j,i)-un(j,i)*(dt/dx)*(vn(j,i)-vn(j,i-1))...
                -vn(j,i)*(dt/dy)*(vn(j,i)-vn(j-1,i))...
                -(dt/(rho*2*dy))*(pn(j+1,i)-pn(j-1,i))...
                +nu*((dt/dx^2)*(vn(j,i+1)-2*vn(j,i)+vn(j,i-1))...
                +(dt/dy^2)*(vn(j+1,i)-2*vn(j,i)+vn(j-1,i)));
            
            b(j,i)=(1/(2*dt))*(((un(j,i+1)-un(j,i-1))/dx)+((vn(j+1,i)-vn(j-1,i))/dy))...
                -((un(j,i+1)-un(j,i-1))/(2*dx))*((un(j,i+1)-un(j,i-1))/(2*dx))...
                -2*((un(j+1,i)-un(j-1,i))/(2*dy))*((vn(j,i+1)-vn(j,i-1))/(2*dx))...
                -((vn(j+1,i)-vn(j-1,i))/(2*dy))*((vn(j+1,i)-vn(j-1,i))/(2*dy));
                                              
            p(j,i)=((pn(j,i+1)+pn(j,i-1))*(dy^2)+(pn(j+1,i)+pn(j-1,i))*(dx^2))...
                /(2*(dx^2+dy^2))-(rho*(dx^2)*(dy^2)*b(j,i))/(2*(dx^2+dy^2));

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
    end
end

%Output
contourf(x,y,p);
hold on
quiver(x,y,u,v,3);
%streamline(x,y,u,v);
hold off
toc                 
