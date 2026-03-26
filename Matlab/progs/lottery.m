% lottery.m
% function x=lottery(k,n)
% draw k lots from a total of n  (k<=n)
% if k==n: (random) reshuffle of (1:n)
% A. van Oosterom; 20110331

function x=lottery(k,n)
if k>n,'k should not exceed n', pause, end
list=1:n;
draw=ceil(n*rand(1,1));
x=draw;
list(list==draw)=[];
n=n-1;
while size(x,2)<k,
      draw=ceil(n*rand(1,1));
      x=[x list(draw)];
      list(draw)=[];
      n=n-1;
end
    