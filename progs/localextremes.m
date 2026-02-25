% local extremes.m
% EXTR=localextremes(ITRI,fun)
% finds local extremes EXTR of a scalar function fun defined at the nodes VER of
% a triangular mesh
% lists: for each local extreme:
% for local maximum:     EXTR(:,[vertex value  1])
% for local flat plane:  EXTR(:,[vertex value  0])
%           minimum:     EXTR(:,[vertex value -1]) 

% A. van Oosterom; 2005/04/06; VER no longer involved  

function EXTR=localextremes(ITRI,fun)
nvals=length(fun);
EXTR=[];
for node=1:nvals,
    bver=buren(ITRI,node);
    bver=bver(bver~=node);
    if       sum(fun(node)<=fun(bver))==0, EXTR=[EXTR; node fun(node)  1];,
      elseif sum(fun(node)~=fun(bver))==0, EXTR=[EXTR; node fun(node)  0];   
      elseif sum(fun(node)>=fun(bver))==0, EXTR=[EXTR; node fun(node)  -1];
    end
end

