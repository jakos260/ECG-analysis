% findrot.m

% interactive selection of desired rotation of at triangulated surface
% should be preceded by a call to triplot
% 20040304; A. van Oosterom


if ~exist('shift'), shift=[0 0 0]; end
phi=0;
theta=0;
gamma=0;

slphibox=[0 .3 .3 .03];
slphi=uicontrol('style','slider','units','norm','position',slphibox,...
'min',-180,'max',180,'Value',0,'CallBack','setphi');

slthetabox=[0 .25 .3 .03];
sltheta=uicontrol('style','slider','units','norm','position',slthetabox,...
'min',0,'max',180,'Value',0,'CallBack','settheta');

slgammabox=[0 .2 .3 .03];
slgamma=uicontrol('style','slider','units','norm','position',slgammabox,...
'min',-180,'max',180,'Value',0,'CallBack','setgamma');


boxt1=[.3 .29 .06 .04];
uit1=uicontrol('style','text','units','norm','position',boxt1);
set(uit1,'string',num2str(phi));

boxt2=[.3 .24 .06 .04];
uit2=uicontrol('style','text','units','norm','position',boxt2);
set(uit2,'string',num2str(theta));

boxt3=[.3 .19 .06 .04];
uit3=uicontrol('style','text','units','norm','position',boxt3);
set(uit3,'string',num2str(gamma));

