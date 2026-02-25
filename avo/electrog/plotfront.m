% plotfront.m
% script of electrog
% plot planar wave front

clf
xb=b*sin(2*phi);
yb=b*cos(2*phi);
xa=a*sin(2*phi);
ya=a*cos(2*phi);
plot(xb-.5,yb)
axis([-1.5 1.5 -1.5 1.5])
axis square
hold on
plot(xa-.5,ya,'r')

for i=2:ntims-1,
  y1=b*cos(phib(i));
  x1=b*sin(phib(i));
  y2=a*cos(phia(i));
  x2=a*sin(phia(i));
  if y1==b, y1=a+x2;end
  if y2==a, y2=b-x1;end
  if y2==-a, y2=min([-b+x1 -a]);end
  plot([x1 x2]-.5,[y1 y2],'r')
  plot([-x1 -x2]-.5,[y1 y2],'r')
end

axis square
if stim=='epi ', title('epicardial stimulus'),end
if stim=='endo', title('endocardial stimulus'), end

setuis