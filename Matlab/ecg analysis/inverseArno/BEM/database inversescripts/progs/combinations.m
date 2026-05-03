% combinations.m
% function ncom=combinations(n,k)
% find number of k combinations taken from n elements
% 070128
function ncom=combinations(n,k)
ncom=1;
if n>1 & k>-1 & k<=n,
   ncom=1;
   if k>n/2, k=n-k; end
   if k>0,
       
       for i=1:k,
          ncom=ncom*(n-i+1)/i;
       end
   end
else,
    'paused because of invalid input'
    pause
end
