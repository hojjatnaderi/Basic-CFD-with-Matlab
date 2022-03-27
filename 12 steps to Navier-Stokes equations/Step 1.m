clc;
clear;
nx=41;
dx=2/(nx-1);
nt=25;
dt=0.025;
c=1;
x=(0:dx:2);
u=ones(1,nx);%u_zero
for i=1:nx
    if (0.5<=x(i)) && (x(i)<=1)
        u(i)=2;
    end
end

un=ones(1,nx);%un is middle parameter and is "u in time"

for n=1:nt
    un=u;%setting u_zero to un
    for i=2:nx
    u(i)=un(i)-c*(dt/dx)*(un(i)-un(i-1));
    end
end

