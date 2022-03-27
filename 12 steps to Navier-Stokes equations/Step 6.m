%Step 6: 2-D Convection (using "array operations" instead of "loops")
%Look at the end of this code for inefficient "loops" method
clc
close all
clear
tic

%Problem characterstics 
nx=101;
ny=101;
nt=80;
sigma=0.2;
dx=(2-0)/(nx-1);
dy=(2-0)/(ny-1);
dt=dx*sigma;

%Initial conditions
u=ones(ny,nx);
u(int64(.5/dy):int64(1/dy+1),int64(.5/dx):int64(1/dx+1))=2;
v=ones(ny,nx);
v(int64(.5/dy):int64(1/dy+1),int64(.5/dx):int64(1/dx+1))=2;
uzero=u;
vzero=v;

%Go-between parameters
un=ones(ny,nx);
vn=ones(ny,nx);

%Iteration
for n=1:nt
    un=u;
    vn=v;
    u(2:end,2:end)=un(2:end,2:end)-un(2:end,2:end)*(dt/dx).*(un(2:end,2:end)-un(2:end,1:end-1))... %term .* is absoloutely important
        -vn(2:end,2:end)*(dt/dy).*(un(2:end,2:end)-un(1:end-1,2:end));
    
    v(2:end,2:end)=vn(2:end,2:end)-un(2:end,2:end)*(dt/dx).*(vn(2:end,2:end)-vn(2:end,1:end-1))...
        -vn(2:end,2:end)*(dt/dy).*(vn(2:end,2:end)-vn(1:end-1,2:end));
    
    %below equations aren't needed becasue of u,v=ones(ny,nx);
    %u(1,:)=1;
    %u(end,:)=1;
    %u(:,1)=1;
    %u(:,end)=1;
    
    %v(1,:)=1;
    %v(end,:)=1;
    %v(:,1)=1;
    %v(:,end)=1;
        
end

%Output illustration
x=(0:dx:2);
y=(0:dy:2);
surf(x,y,u);
toc

%Below is inefficient "loops" method. 
%
%Iteration
%for n=1:nt
%    un=u;
%    vn=v;
%    for j=2:ny
%        
%        for i=2:nx
%            u(j,i)=un(j,i)-un(j,i)*(dt/dx)*(un(j,i)-un(j,i-1))...
%                -vn(j,i)*(dt/dy)*(un(j,i)-un(j-1,i));
%            
%            v(j,i)=vn(j,i)-un(j,i)*(dt/dx)*(vn(j,i)-vn(j,i-1))...
%                -vn(j,i)*(dt/dy)*(vn(j,i)-vn(j-1,i));
%            
%             u(1,:)=1;
%             u(end,:)=1;
%             u(:,1)=1;
%             u(:,end)=1;
%    
%             v(1,:)=1;
%             v(end,:)=1;
%             v(:,1)=1;
%             v(:,end)=1;
%        end
%    end
%        
%end
