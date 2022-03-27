function [pout]=laplace2d(p,y,dx,dy,l1norm_target)
%this function solves 2D Laplace equation

l1norm=1;
while l1norm>l1norm_target
    pn=p;
    p(2:end-1,2:end-1)=((dy^2)*(pn(2:end-1,3:end)+pn(2:end-1,1:end-2))+...
        (dx^2)*(pn(3:end,2:end-1)+p(1:end-2,2:end-1)))...
        /(2*(dx^2+dy^2));
    
    %Boundary conditions
    p(:,1)=0;
    p(:,end)=y;
    p(1,:)=p(2,:);
    p(end,:)=p(end-1,:);

    %Ee want to sum all elements in a matrix
    l1norm=sum(sum((abs(p)-abs(pn))))/sum(sum(abs(pn))); 
end
pout=p;
end
