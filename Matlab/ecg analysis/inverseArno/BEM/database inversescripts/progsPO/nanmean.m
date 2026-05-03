function m=nanmean(x,dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Name:    nanmean
%
% Desc:    calculate mean ignoring nan
%
% Inputs:  X, dim see mean
%
% Outputs: m: mean
%
% Created: 13-Nov-2013 08:45:54
%
% Author:  peteroosterhoff
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inan=isnan(x);
x(inan)=0;
if ~exist('dim','var')
    n=sum(~inan);
    n(n==0)=NaN;
    m=sum(x)./n;
else
    n=sum(~inan,dim);
    n(n==0)=NaN;
    m=sum(x,dim)./n;
end
