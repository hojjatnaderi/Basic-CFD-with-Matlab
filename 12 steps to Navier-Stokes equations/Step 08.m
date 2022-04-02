%Step 8: Burgers' Equation in 2D
%This method uses faster "array operations" technique.
%Down below (after this code and in comments) I'll put the slower "interloops" method. 
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
        
    u(2:end-1,2:end-1)=un(2:end-1,2:end-1)...
        -(dt/dx)*un(2:end-1,2:end-1).*(un(2:end-1,2:end-1)-un(2:end-1,1:end-2))... %remove dot before * and see what happens
        -(dt/dy)*vn(2:end-1,2:end-1).*(un(2:end-1,2:end-1)-un(1:end-2,2:end-1))...
        +nu*(dt/(dx)^2)*(un(2:end-1,3:end)-2*un(2:end-1,2:end-1)+un(2:end-1,1:end-2))...
        +nu*(dt/(dy)^2)*(un(3:end,2:end-1)-2*un(2:end-1,2:end-1)+un(1:end-2,2:end-1));
    
    v(2:end-1,2:end-1)=vn(2:end-1,2:end-1)...
        -(dt/dx)*un(2:end-1,2:end-1).*(vn(2:end-1,2:end-1)-vn(2:end-1,1:end-2))...
        -(dt/dy)*vn(2:end-1,2:end-1).*(vn(2:end-1,2:end-1)-vn(1:end-2,2:end-1))...
        +nu*(dt/(dx)^2)*(vn(2:end-1,3:end)-2*vn(2:end-1,2:end-1)+vn(2:end-1,1:end-2))...
        +nu*(dt/(dy)^2)*(vn(3:end,2:end-1)-2*vn(2:end-1,2:end-1)+vn(1:end-2,2:end-1));
    
    %Boundary Conditions:
    u(1,:)=1;
    u(end,:)=1;
    u(:,1)=1;
    u(:,end)=1;
    
    v(1,:)=1;
    v(end,:)=1;
    v(:,1)=1;
    v(:,end)=1;

end
%
x=(0:dx:2);
y=(0:dy:2);
surf(x,y,u);
toc
% For comparing purposes, below I'll put the inefficient code
% which uses "interloops" instead of "array operations". 
%
%Slower Method:
%2d Burgers equation
%clc
%close all
%clear
%tic
%
%nx=41;
%ny=41;
%nt=120;
%dx=2/(nx-1);
%dy=2/(ny-1);
%sigma=0.0009;
%nu=0.01;
%dt=sigma*dx*(dy/nu);
%
%u=ones(ny,nx);
%v=ones(ny,nx);
%u(int64(.5/dy):int64(1/dy+1),int64(.5/dx):int64(1/dx+1))=2;
%v(int64(.5/dy):int64(1/dy+1),int64(.5/dx):int64(1/dx+1))=2;
%
%for n=1:nt
%    un=u;
%    vn=v;
%    for j=2:ny-1
%       for i=2:nx-1
%            u(j,i)=un(j,i)-(dt/dx)*un(j,i)*(un(j,i)-un(j,i-1))...
%                -(dt/dy)*vn(j,i)*(un(j,i)-un(j-1,i))...
%                +nu*(dt/(dx)^2)*(un(j,i+1)-2*un(j,i)+un(j,i-1))...
%               +nu*(dt/(dy)^2)*(un(j+1,i)-2*un(j,i)+un(j-1,i));
%            
%            v(j,i)=vn(j,i)-(dt/dx)*un(j,i)*(vn(j,i)-vn(j,i-1))...
%                -(dt/dy)*vn(j,i)*(vn(j,i)-vn(j-1,i))...
%                +nu*(dt/(dx)^2)*(vn(j,i+1)-2*vn(j,i)+vn(j,i-1))...
%                +nu*(dt/(dy)^2)*(vn(j+1,i)-2*vn(j,i)+vn(j-1,i));
%            
%            u(1,:)=1;
%            u(end,:)=1;
%            u(:,1)=1;
%            u(:,end)=1;
%            
%            v(1,:)=1;
%            v(end,:)=1;
%            v(:,1)=1;
%            v(:,end)=1;
%        end
%    end
%end
%
%x=(0:dx:2);
%y=(0:dy:2);
%surf(x,y,u);
%toc            
