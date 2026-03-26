% rflin0.m
% straight line with zero intercept
% function [Yest,G]=rflin1(X,p,mode);
% y=p(0)*x

% funtype=1;


    function [Yest,G]=rflin1(X,p,mode);
    [nx ndum]=size(X);
    [np ndum]=size(p);
    Yest=p(1)*X;
    
    if mode==0, return, end
    G(1:nx,1)=X;
    