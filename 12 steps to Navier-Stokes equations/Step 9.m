%2D Laplace equation
clc
close all
clear
tic %time measurement

%Setting up the field
nx=31;
ny=31;
dx=2/(nx-1);
dy=1/(ny-1);
x=(0:dx:2);
y=(0:dy:1);
l1norm_target=1e-4;

%Initial conditions
p=zeros(ny,nx);

%Convergence parameters
l1norm=1;
n=1;

%Iteration
while l1norm>l1norm_target
    pn=p;
    p(2:end-1,2:end-1)=((dy^2)*(pn(2:end-1,3:end)+pn(2:end-1,1:end-2))+...
        (dx^2)*(pn(3:end,2:end-1)+pn(1:end-2,2:end-1)))...
        /(2*(dx^2+dy^2));

    p(:,1)=0;
    p(:,end)=y;
    p(1,:)=p(2,:);
    p(end,:)=p(end-1,:);

    %we want to sum all elements in a matrix
    l1norm=sum(sum((abs(p)-abs(pn))))/sum(sum(abs(pn))); 
    n=n+1;
end

%Output
surf(x,y,p)
toc
