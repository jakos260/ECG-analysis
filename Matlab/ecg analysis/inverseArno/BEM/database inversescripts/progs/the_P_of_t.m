function p = the_P_of_t(t,df);
%function p = the_P_of_t(t,df);
% rechter overschrijdingskans van t_distribution met df degrees of freedom
% A. van Oosterom; 2007_12_13

t=abs(t);
df=ceil(abs(df));

% f=1/(sqrt(df)*B(0.5,df/2))* int_0^t  (1+x^2/df)^(-(df+1)/2) )
% B(0.5,df/2)= (complete) 
% Betafunction= Gamma(0.5)*Gamma(df/2)/(Gamma((df+1)/2)
% Abramowitz;  page 948
 
nx=401;
fnorm=sqrt(df*pi)*gamma_function(df/2)/gamma_function((df+1)/2);
i=0:nx-1;
x=i*t/(nx-1);
f=(1+x.^2/df).^(-(df+1)/2)/fnorm;

% figure(1)
% clf
% plot(x,f)
p=0.5-max(introws(f))*t/(nx-1);


