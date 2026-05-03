% file ttest.m
% [t,df]=ttest(x,y,mode)
% compute t value of various hypotheses
% default: mode==0: t_test, testing mu_x=y; (default y=0)  sigma_x unknown
%          mode==1:       , testing mu_x==mu_y; sigmas equal but unknown
%          mode==2;       , paired t_test: testing mu_x==mu_y; sigmas equal but unknown
%          mode==3;       , t-statistic testing zero correlation between x and y  

function [t,df]=ttest(x,y,mode)
nx=length(x);
t=0; df=0;

if nargin==1, 
  mux=0;
end

if nargin==2, 
   mux=y;
end

if nargin<=2,
    mx=mean(x-mux);
    sx=std(x);
    t=mx*sqrt(nx)/sx;
    df=nx-1;
end

if nargin>=2 & mode >0,
   if mode==1,
       nx=length(x);
       ny=length(y);
       sx=std(x);
       sy=std(y);
       s=sqrt(((nx-1)*sx^2+(ny-1)*sy^2)/(nx+ny-2));
       t=(mean(x)-mean(y))*sqrt(nx*ny)/(s*sqrt(nx+ny));
       df=nx+ny-2;
   end
   if mode==2,
       % paired t_test
       if length(x)~=length(y), 'error: x and y have unequal length', end
       t=mean(x-y)*sqrt(nx)/std(x-y);
       df=nx-1;
   end
   if mode==3,
       % t_test; testing zero correlation between x and y
       % Basic&Clinical Biostatistics; B.Dawson-Saunders and R. G. Trapp;
       % Appleton and Lange; Norwalk, Connecticut; 1994; page 165
       if length(x)~=length(y), 'error: x and y have unequal length', end
       RHO=corrcoef(x,y);
       r=RHO(1,2);
       t=r*sqrt(nx-2)/sqrt(1-r^2);
       df=nx-2;
   end
end
   
       
       



