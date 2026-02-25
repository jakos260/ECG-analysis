% function i=icyc(k,n)
% k and n are integers;  n a scalar
% 20110204
% rem type of function; example:
% for n=3 and 
% k=[ -3 -2 -1 0 1 2 3 4 5 6 7 etc ]; it produces
% i=[  3  1  2 3 1 2 3 1 2 3 1 etc;;
function k=icyc(k,n)
  k=rem(k,n);
  k(k<=0)= k(k<=0)+n;
end