
% function gam=gamma_function(z)
% n=floor(z); x=z-n;
% n integer; (n> 0;  x=0)   or   (n>=0; x=0.5);
% for x=0, gam=(n-1)!% for x=0; the function returns (n-1)!
% for more values of x: see Abramowitz page 255

% A. van Oosterom; 2007_12_13

function gam=gamma_function(z);

n=floor(z); x=z-n;
n=floor(abs(n));
if x~=0; x=0.5; end
if n+x < 0.5, 'non-supported input values', pause, end 

if  x==0,
    gam=1;
else,
    gam=sqrt(pi);
end

if  n+x < 0.5, return, end

if x==0,
   y=1;
   for i=2:n,
       gam=y*gam;
       y=y+1;
   end
else,
    y=x-1; 
    for i=1:n,
        y=y+1;
        gam=y*gam;
        
    end
end








