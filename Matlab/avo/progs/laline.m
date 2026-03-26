% laline.m
% [lambda, dis]=laline(r,r1,r2)
% lambda: specifies projection of vector r on a line segment
% from r1 to r2
% dis=distance from r to line segment

function [lambda, dis]=laline(r,r1,r2)
  s1=r2-r1;
  s2= r-r1;
  s1dot=s1*s1';
  lambda=0;
  dis=0;
  s2dot=s2*s2';

if s2dot > eps,
  s2dot=sqrt(s2dot);
  sdot=s1*s2';
  lambda=sdot/s1dot;
  s1dot=sqrt(s1dot);
  cosine=sdot/(s1dot*s2dot);  
  d=sqrt(1-cosine^2)*s2dot;
   dis=d;
end

end     