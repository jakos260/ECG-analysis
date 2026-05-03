% example 
clear
clf
p = randn(1000,1);
x=[-4:.4:4];
[N, X] = hist(p,x);
h = bar(X, N, 'hist');
set(h,'Facecolor','none','Edgecolor','r')
hold on
p = randn(1000,1);
[N, X] = hist(p,x);
h = bar(X, N, 'hist');
set(h,'Facecolor','none','Edgecolor','b')
for i=1:101,
x(i)=-4+8*(i-1)/100;
y(i)=400/sqrt(2*pi)*exp(-x(i)^2/2);
end
plot(x,y,'k')



%ch = get(gca,'ch');
%h = bar(X, N, 'histc');
%    _____
%    )o|o(
%-oO0--V--0Oo--------------------------------------------
%Ali Taylan Cemgil        http://www.mbfys.kun.nl/~cemgil 
%SNN    -    University of Nijmegen,      The Netherlands

