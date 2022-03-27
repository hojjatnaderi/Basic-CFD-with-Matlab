%Navier Stokes - Channel flow - HJ - Last step
%this code doesn't iterate pn or p
clc
close all
clear
tic %time measurement

%Setting up the field 
nx=41;
ny=41;
dx=2/(nx-1);
dy=2/(ny-1);

x=(0:dx:2);
y=(0:dy:2);

rho=1;
nu=0.1;
f=1;
dt=0.01;

%Initial conditions
u=zeros(ny,nx);
v=zeros(ny,nx);
p=zeros(ny,nx);
b=zeros(ny,nx);

un=zeros(ny,nx);
vn=zeros(ny,nx);
pn=zeros(ny,nx);

%Convergence parameters
udiff=1;
stepcount=0;

%Iteration
while udiff>0.00001
    un=u;
    vn=v;
    pn=p;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pressure equation
    
        b(2:end-1,2:end-1)=(1/(2*dt))*(((un(2:end-1,3:end)-un(2:end-1,1:end-2))/dx)+((vn(3:end,2:end-1)-vn(1:end-2,2:end-1))/dy))...
            -((un(2:end-1,3:end)-un(2:end-1,1:end-2))/(2*dx)).*((un(2:end-1,3:end)-un(2:end-1,1:end-2))/(2*dx))...
            -2*((un(3:end,2:end-1)-un(1:end-2,2:end-1))/(2*dy)).*((vn(2:end-1,3:end)-vn(2:end-1,1:end-2))/(2*dx))...
            -((vn(3:end,2:end-1)-vn(1:end-2,2:end-1))/(2*dy)).*((vn(3:end,2:end-1)-vn(1:end-2,2:end-1))/(2*dy));

        p(2:end-1,2:end-1)=((pn(2:end-1,3:end)+pn(2:end-1,1:end-2))*(dy^2)+(pn(3:end,2:end-1)+pn(1:end-2,2:end-1))*(dx^2))...
            /(2*(dx^2+dy^2))-(rho*(dx^2)*(dy^2)*b(2:end-1,2:end-1))/(2*(dx^2+dy^2));

        %periodic boundary condition (B.C.) at x=2
        b(2:end-1,end)=(1/(2*dt))*(((un(2:end-1,1)-un(2:end-1,end-1))/dx)+((vn(3:end,end)-vn(1:end-2,end))/dy))...
            -((un(2:end-1,1)-un(2:end-1,end-1))/(2*dx)).*((un(2:end-1,1)-un(2:end-1,end-1))/(2*dx))...
            -2*((un(3:end,end)-un(1:end-2,end))/(2*dy)).*((vn(2:end-1,1)-vn(2:end-1,end-1))/(2*dx))...
            -((vn(3:end,end)-vn(1:end-2,end))/(2*dy)).*((vn(3:end,end)-vn(1:end-2,end))/(2*dy));

        %periodic B.C. at x=0
        b(2:end-1,1)=(1/(2*dt))*(((un(2:end-1,2)-un(2:end-1,end))/dx)+((vn(3:end,1)-vn(1:end-2,1))/dy))...
            -((un(2:end-1,2)-un(2:end-1,end))/(2*dx)).*((un(2:end-1,2)-un(2:end-1,end))/(2*dx))...
            -2*((un(3:end,1)-un(1:end-2,1))/(2*dy)).*((vn(2:end-1,2)-vn(2:end-1,end))/(2*dx))...
            -((vn(3:end,1)-vn(1:end-2,1))/(2*dy)).*((vn(3:end,1)-vn(1:end-2,1))/(2*dy));


        %periodic B.C. at x=2
        p(2:end-1,end)=((pn(2:end-1,1)+pn(2:end-1,end-1))*(dy^2)+(pn(3:end,end)+pn(1:end-2,end))*(dx^2))...
            /(2*(dx^2+dy^2))-(rho*(dx^2)*(dy^2)*b(2:end-1,end))/(2*(dx^2+dy^2));

        %periodic B.C. at x=0
        p(2:end-1,1)=((pn(2:end-1,2)+pn(2:end-1,end))*(dy^2)+(pn(3:end,1)+pn(1:end-2,1))*(dx^2))...
            /(2*(dx^2+dy^2))-(rho*(dx^2)*(dy^2)*b(2:end-1,1))/(2*(dx^2+dy^2));

        %wall B.C. at y=0 & y=2 we have dp/dy=0
        p(1,:)=p(2,:);
        p(end,:)=p(end-1,:);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% velocities equations
        %                                             Using a dot before * is optional (.* is optional)
        u(2:end-1,2:end-1)=un(2:end-1,2:end-1)...
            -un(2:end-1,2:end-1)*(dt/dx).*(un(2:end-1,2:end-1)-un(2:end-1,1:end-2))...
            -vn(2:end-1,2:end-1)*(dt/dy).*(un(2:end-1,2:end-1)-un(1:end-2,2:end-1))...
            -(dt/(rho*2*dx))*(pn(2:end-1,3:end)-pn(2:end-1,1:end-2))...
            +nu*((dt/dx^2)*(un(2:end-1,3:end)-2*un(2:end-1,2:end-1)+un(2:end-1,1:end-2))...
            +(dt/dy^2)*(un(3:end,2:end-1)-2*un(2:end-1,2:end-1)+un(1:end-2,2:end-1)))+dt*f;
            
        v(2:end-1,2:end-1)=vn(2:end-1,2:end-1)...
            -un(2:end-1,2:end-1)*(dt/dx).*(vn(2:end-1,2:end-1)-vn(2:end-1,1:end-2))...
            -vn(2:end-1,2:end-1)*(dt/dy).*(vn(2:end-1,2:end-1)-vn(1:end-2,2:end-1))...
            -(dt/(rho*2*dy))*(pn(3:end,2:end-1)-pn(1:end-2,2:end-1))...
            +nu*((dt/dx^2)*(vn(2:end-1,3:end)-2*vn(2:end-1,2:end-1)+vn(2:end-1,1:end-2))...
            +(dt/dy^2)*(vn(3:end,2:end-1)-2*vn(2:end-1,2:end-1)+vn(1:end-2,2:end-1)));
            
        %periodic B.C. x=2                             .* is essential
        u(2:end-1,end)=un(2:end-1,end)-un(2:end-1,end)*(dt/dx).*(un(2:end-1,end)-un(2:end-1,end-1))...
            -vn(2:end-1,end)*(dt/dy).*(un(2:end-1,end)-un(1:end-2,end))...   
            -(dt/(rho*2*dx))*(pn(2:end-1,1)-pn(2:end-1,end-1))...
            +nu*((dt/dx^2)*(un(2:end-1,1)-2*un(2:end-1,end)+un(2:end-1,end-1))...
            +(dt/dy^2)*(un(3:end,end)-2*un(2:end-1,end)+un(1:end-2,end)))+dt*f;

        %periodic B.C. x=0
        u(2:end-1,1)=un(2:end-1,1)-un(2:end-1,1)*(dt/dx).*(un(2:end-1,1)-un(2:end-1,end))...
            -vn(2:end-1,1)*(dt/dy).*(un(2:end-1,1)-un(1:end-2,1))...
            -(dt/(rho*2*dx))*(pn(2:end-1,2)-pn(2:end-1,end))...
            +nu*((dt/dx^2)*(un(2:end-1,2)-2*un(2:end-1,1)+un(2:end-1,end))...
            +(dt/dy^2)*(un(3:end,1)-2*un(2:end-1,1)+un(1:end-2,1)))+dt*f;
            
        %periodic B.C. x=2 for v
         v(2:end-1,end)=vn(2:end-1,end)-un(2:end-1,end)*(dt/dx).*(vn(2:end-1,end)-vn(2:end-1,end-1))...
            -vn(2:end-1,end)*(dt/dy).*(vn(2:end-1,end)-vn(1:end-2,end))...
            -(dt/(rho*2*dy))*(pn(3:end,end)-pn(1:end-2,end))...
            +nu*((dt/dx^2)*(vn(2:end-1,1)-2*vn(2:end-1,end)+vn(2:end-1,end-1))...
            +(dt/dy^2)*(vn(3:end,end)-2*vn(2:end-1,end)+vn(1:end-2,end)));
            
        %periodic B.C. x=0 for v
        v(2:end-1,1)=vn(2:end-1,1)-un(2:end-1,1)*(dt/dx).*(vn(2:end-1,1)-vn(2:end-1,end))...
            -vn(2:end-1,1)*(dt/dy).*(vn(2:end-1,1)-vn(1:end-2,1))...
            -(dt/(rho*2*dy))*(pn(3:end,1)-pn(1:end-2,1))...
            +nu*((dt/dx^2)*(vn(2:end-1,2)-2*vn(2:end-1,1)+vn(2:end-1,end))...
            +(dt/dy^2)*(vn(3:end,1)-2*vn(2:end-1,1)+vn(1:end-2,1)));
            
        %bottom and top wall B.C.
        u(1,:)=0;
        u(end,:)=0;

        v(1,:)=0;
        v(end,:)=0;

    udiff=(sum(u)-sum(un))/sum(u);
    stepcount=stepcount+1;
end

%Making a third of output visible
xsols=x(1:3:end);
ysols=y(1:3:end);
usols=u(1:3:end,1:3:end);
vsols=v(1:3:end,1:3:end);

%Output illustration
%contourf(x,y,u);
hold on
quiver(xsols,ysols,usols,vsols,1);
%streamline(x,y,u,v);
hold off
toc            
