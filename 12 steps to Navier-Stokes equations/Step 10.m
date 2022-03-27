%Poisson's equation
clc
close all
clear
tic %time measurement

%Setting up the field
nx=51;
ny=51;
nt=100;
xmin=0;
xmax=2;
ymin=0;
ymax=1;

dx=(xmax-xmin)/(nx-1);
dy=(ymax-ymin)/(ny-1);

%Initial conditions
p=zeros(ny,nx);
pd=zeros(ny,nx);
b=zeros(ny,nx);
b(int64(ny/4),int64(nx/4))=100;
b(int64(3*ny/4),int64(3*nx/4))=-100;

%Iteration
for it=1:nt
    pd=p;
    p(2:end-1,2:end-1)=((dy^2)*(pd(2:end-1,3:end)+pd(2:end-1,1:end-2))+...
        (dx^2)*(pd(3:end,2:end-1)+pd(1:end-2,2:end-1))...
        -b(2:end-1,2:end-1)*(dx^2)*(dy^2))...
        /(2*(dx^2+dy^2));
        
        %Boundary conditions
        p(1,:)=0;
        p(end,:)=0;
        p(:,1)=0;
        p(:,end)=0;        
end

%Output illustration
x=(xmin:dx:xmax);
y=(ymin:dy:ymax);
surf(x,y,p);
toc  
