% zoomin.m
% script of select time segment
tbegold=tbeg;
interval=tend-tbeg;
tbeg=tbeg+round(interval/4);
tb=max(tb,tbeg);


tend=tend-round(interval/4);
te=min(te,tend);

column=max(column+tbegold-tbeg,1); 
column=min(column,te);
% 
% set(text5,'string',num2str(tb))
% set(text3,'string',num2str(te))

setplots
