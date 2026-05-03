% rflogist.m
% function [y,G]=rflogist(x,p,mode);
% y=p(1)/(1+exp(-p(2)*(x-p(3)))

% funtype=3

    function [y,G]=rflogist(x,p,mode);
    [nx ndum]=size(x);
    [np ndum]=size(p);
    ecx=exp(-p(2)*(x-p(3)));
    tmp=1./(1+ecx);
    y=p(1)*tmp; 
    
    %if mode==0, return, end
    G(:,1)=tmp;
    tmp=tmp.^2;
    G(:,2)= p(1)*tmp.*ecx.*(x-p(3));
    G(:,3)=-p(1)*p(2)*tmp.*ecx*p(2);

