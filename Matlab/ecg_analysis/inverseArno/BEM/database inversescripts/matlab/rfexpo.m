% rfexpo.m
% function [y,G]=rfexpo(x,p,mode);
% y=p(1)+exp(p(2)*(x-p(3))

% funtype=4

    function [y,G]=rfexpo(x,p,mode);
    [nx ndum]=size(x);
    [np ndum]=size(p);
    ecx=exp(p(2)*(x-p(3)));
    y=p(1)+ecx; 
    
    if mode==0, return, end
    
    G(1:nx,1)=1;
    G(:,2)=ecx.*(x-p(3));
    G(:,3)=-ecx*p(2);


