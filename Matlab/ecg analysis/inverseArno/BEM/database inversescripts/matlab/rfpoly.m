% rfpoly.m
% function [y,G]=rfpoly(x,p,mode);
% y=polynomial of x; G the derivatives  
% y and x column vectors
% y=p(1)+p(2)*x+p(3)*x^2 etc

% funtype= 2

    function [y,G]=rfpoly(x,p,mode);
    [ni jdum]=size(x);
    np=length(p);
    k=np;
    Y=[];
    % apply Horners rule

    if np>1,
        y=p(np)*x;
        k=k-1;
        while k>1,
           y=x.*(p(k)+y);
           k=k-1;
         end
     end
     y=y+p(1);
    
    if mode==0, return, end
    G=zeros(ni,np);
    xpower=ones(ni,1);
    G(:,1)=xpower;
    for k=2:np;
        xpower=xpower.*x;
        G(:,k)=xpower;
    end
    