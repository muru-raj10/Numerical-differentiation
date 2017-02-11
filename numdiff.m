f=@(x) exp(x);
h=0.01;
x0=0;
v=[-1/2 1 0 -1 1/2];
f3 = sum(v.*f(x0+(-2:2)*h))/h^3;

d=[-11/6 3 -3/2 1/3];
f1=sum(d.*f(x0+(0:3)*h))/h;

%this gives the central diff table
format rat
A=[1 1 1 1 1;
(-2:2);
(-2:2).^2/factorial(2);
(-2:2).^3/factorial(3);
(-2:2).^4/factorial(4)];

b=[0;1;0;0;0];
inv(A)*b


