function linconv(nx)
%this function solves the linear convection equation 
%Setting up the field
dx=(2-0)/(nx-1);
nt=20;
c=1;
sigma=0.5; %Courant number
dt=sigma*dx; %instead of separate definition of dt & dx
%Defining u_zero
u=ones(1,nx);
u(int64((0.5-0)/dx): int64((1-0)/dx+1))=2;
%Plotting u_zero
u_zero=u;
x=(0:dx:2);
%plot(x,u_zero)
%hold on
%Defining un parameter
un=ones(1,nx);
%Iterating
 for n=1:nt
     un=u;
     for i=2:nx
         u(i)=un(i)-c*(dt/dx)*(un(i)-un(i-1));
     end
 end
 %Output
 plot(x,u)
end
