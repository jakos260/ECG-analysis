% minpar.m
% determine x-value for which f(x) is minimal on the interval
% x(1) < x(3)
% local parabole assumed; function specified at
% x(1) < x(2) < x(3)
% specified by f(x(1)), f(x(2)) and f(x(3))
% because of its application in an optimization procedure (finding the 
% optimal step size; a positive value) the value xmin is put at x(3)
% in case of a horizontal line
function xmin=minpar(f,x)
x=x';
f=f';
frac=(x(2)-x(1))/(x(3)-x(1));
f2=f(1)*(1-frac)+f(3)*frac;
if f(2) >= f2,
   if f(1)<f(3),
   xmin=x(1);
   else,
   xmin=x(3);
   end
else,
A=[x.^2 x ones(3,1)];
a=pinv(A)*f;
xmin=-a(2)/(2*a(1));
end
