% rf_2expos.m
% function [y,G]=rf_2expos(x,p,mode);
% y=p(1)*exp(p(2)*x) + p(3)*exp(p(4)*(x-max(x))) + p(5)

% funtype=8

    function [y,G]=rf_2expos(x,p,mode);
    [nx ndum]=size(x);
    [np ndum]=size(p);
    
    exp1=exp(p(2)*x);
    exp2=exp(p(4)*(x-max(x)));
    y=p(1)*exp1+p(3)*exp2+p(5); 
    
    if mode==0, return, end
    
    G(1:nx,1)=exp1;
    G(:,2)=p(1)*exp1.*x;
    G(:,3)=exp2;
    G(:,4)=p(3)*exp2.*(x-max(x));
    G(:,5)=ones(np,1);
