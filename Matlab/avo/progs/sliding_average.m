function [xsla,ysla]=sliding_average(x,y,winfrac)
% function [xsla,ysla]=sliding_average(x,y,winfrac)
% computes average over slinding x-windo with width: win_frac*range of x
nx=length(x);
if size(x,1)<size(x,2), x=x'; end
if size(y,1)<size(y,2), y=y'; end

xmax=max(x); xmin=min(x);
win=winfrac*(xmax-xmin); winh=win/2;
DATA=[x,y];
DATA=sortrows(DATA,1);
x=DATA(:,1); y=DATA(:,2);
for i=1:nx;
    data=y(x<=x(i)+winh & x>=x(i)-winh);
    %x(i)
    ysla(i)=mean(data);
end
xsla=x;