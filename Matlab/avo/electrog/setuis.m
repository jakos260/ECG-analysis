% setuis.m
% 200301729
% sets all uis

uie1=uicontrol('style','edit');

ui1p1=uicontrol('style','pushbutton');


ui1r1=uicontrol('style','radio');
ui1r2=uicontrol('style','radio');
ui1r3=uicontrol('style','radio');
%ui1r4=uicontrol('style','radio');
%ui1r5=uicontrol('style','radio');
ui1r6=uicontrol('style','radio');
ui1r7=uicontrol('style','radio');
ui1r8=uicontrol('style','radio');
ui1r9=uicontrol('style','radio');

uis1=uicontrol('style','slider');
uis2=uicontrol('style','slider');

uit1=uicontrol('style','text');
uit2=uicontrol('style','text');
uit3=uicontrol('style','text');
uit4=uicontrol('style','text');
uit5=uicontrol('style','text');
uit6=uicontrol('style','text');

set(uie1,'units','norm','position',uiboxe1','string',num2str(a));

set(ui1p1,'units','norm','position',uiboxp1,'string','inner-radius','Callback',...
'getradius')

set(ui1r1,'units','norm','position',uibox1r1,'string','endo',...
'value',r1,'foregroundcolor','r','Callback','endostim')

set(ui1r2,'units','norm','position',uibox1r2,'string','epi',...
'value',r2,'foregroundcolor','b','Callback','epistim')

set(ui1r3,'units','norm','position',uibox1r3,'string','inhom',...
'value',inho,'foregroundcolor','r','Callback','inhom')
%set(ui1r4,'units','norm','position',uibox1r4,'string','epi',...
%'value',r4,'foregroundcolor','b','Callback','partepi')
%set(ui1r5,'units','norm','position',uibox1r5,'string','phi',...
%'value',r5,'foregroundcolor','k','Callback','partall')
set(ui1r6,'units','norm','position',uibox1r6,'string','pots',...
'value',r6,'foregroundcolor','k','Callback','showpots')
set(ui1r7,'units','norm','position',uibox1r7,'string','fronts',...
'value',r7,'foregroundcolor','k','Callback','showfronts')
set(ui1r8,'units','norm','position',uibox1r8,'string','endo',...
'value',r8,'foregroundcolor','r','Callback','obsendo')
set(ui1r9,'units','norm','position',uibox1r9,'string','epi',...
'value',r9,'foregroundcolor','b','Callback','obsepi')

set(uis1,'units','norm','position',uiboxs1,...
'min',.2,'max',2,'Value',sl1,'CallBack','scalpots');
set(uis2,'units','norm','position',uiboxs2,...
'min',0,'max',1.5,'Value',robs,'CallBack','setrobs');

set(uit1,'units','norm','position',uiboxt1,'string','stimulus:')
set(uit2,'units','norm','position',uiboxt2,'string','display:')
set(uit3,'units','norm','position',uiboxt3,'string',num2str(robs));
set(uit4,'units','norm','position',uiboxt4,'string','observe:');
set(uit5,'units','norm','position',uiboxt5,'string','scale');
set(uit6,'units','norm','position',uiboxt6,'string','r-obs');


