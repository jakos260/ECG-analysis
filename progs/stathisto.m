% stathist0.m
% determine mean and standard deviation of 
% data persented by means of a histogram
% assuming bins of equal width
 function [m,s]=stathisto(y)
dim=size(y);
nbins=max(dim);
% correct negative values
for i=1:nbins,
if y(i) < 0 , y(i)=0; end
end
ytot=sum(y);
f=y/ytot;
ex=0;
exsq=0;
for x=1:nbins,
ex=ex+x*f(x);
exsq=exsq+x^2*f(x);
end
m=ex;
s=sqrt(abs(exsq-m^2));
	 