% rfexp0.m
% function [y,G]=rfexp0(x,p,mode);
% y=exp(p(1)*(x-p(2))

% funtype=5

    function [y,G]=rfexp0(x,p,mode);
    [nx ndum]=size(x);
    [np ndum]=size(p);
    ecx=exp(p(1)*(x-p(2)));
    y=ecx; 
    
    if mode==0, return, end
    G(:,1)=ecx.*(x-p(2));
    G(:,2)=-ecx*p(1);


