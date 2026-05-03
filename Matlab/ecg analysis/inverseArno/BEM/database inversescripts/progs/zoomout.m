% zoomout.m
% script of select time segment
tbegold=tbeg;
interval=tend-tbeg;
tbeg=tbeg-round(interval/2);
tbeg=max(1,tbeg);
tend=tbeg+2*interval+1;
tend=min(tend,nt);
column=column+tbegold-tbeg;

if tend==nt, tbeg=max(1,nt-2*interval); end
if tend-tbeg > interval,  
%     set(text5,'string',num2str(tb))
%     set(text3,'string',num2str(te))
    setplots
end