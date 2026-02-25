% minparab.ftn
%C determines x-value for which y(x) is minimal
%C on the interval
%C x(1) < x(3)
%C local parabole assumed, specified at
%C x(1) < x(2) < x(3) by y(1)), y(2)) and y(3))
function xmin=minparab(y,x)
%C check input data
	del12=x(1)-x(2);
	del23=x(2)-x(3);
	del31=x(3)-x(1);
	frac=-del12/del31;
	y2=y(1)*(1-frac)+y(3)*frac
	if y2 <= y(2), 
% straight line detected; minimum at the boundary      
	  if y(1) <= y(3), xmin=x(1);end
% NOTE: handling of equal sign arbitrarily puts xmin at x(1)
	  if y(3) >  y(1), xmin=x(3);end
	else 
% xmin within the interval x(1)--X(3)
	sum=del23*y(1)+del31*y(2)+del12*y(3);
	atel=sum;
	sum=del23*x(1)*x(1)+del31*x(2)*x(2)+del12*x(3)*x(3);
	anoe=sum;
	a=atel/anoe;
	b=(y(1)-y(2)-a*(x(1)*x(1)-x(2)*x(2)))/del12;
	if a ~=0,
	xmin=-b/(2*a);
	end
	if a ==0,
	 if b >=0, xmin=x(1);end
	 if b < 0, xmin=x(3);end
	end
	
	if xmin < x(1), xmin=x(1), end
	if xmin > x(3), xmin=x(3), end
end