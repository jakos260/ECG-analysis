function m=mean1d(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name:    mean1d
% 
% Desc:    calculate mean of x, ignoring NaN's
% 
% Inputs:  x: n-dimensianl matrix, possibly containing NaN's
% 
% Outputs: m: mean of all elements of x not NaN.
% 
% Created: 12-Nov-2013 15:08:04
% 
% Author:  peteroosterhoff
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x(:);
m=mean(x(~isnan(x)));