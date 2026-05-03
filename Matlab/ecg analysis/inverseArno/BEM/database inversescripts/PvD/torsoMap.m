function torsoMap(ni,PHI,av,depV,dir ,fn)

global TORSO
global AS

if ni >0
	PHI=PHI-ones(length(TORSO.TVER),1)*PHI(ni,:);
end
infor1=GetAmplInfos(PHI,max(depV)+av );

hf=figure(2031);set(hf,'PaperPositionMode','auto')
set(gcf,'position',[520   100   1211 406  ]);clf; 
colormap(AS.AVO)

CP=infor1(:,1);
CQ=infor1(:,2);
CT=infor1(:,3);
CST=infor1(:,4);
doampl=1;
doQRSArea=0;
doST=0;
%%
if doampl
hf=figure(2032);set(hf,'PaperPositionMode','auto')
set(gcf,'position',[520   100   1011 406  ]);clf; 
colormap(AS.AVO)
axes('Position',[1-1/4*(4)+0.0 0.01 0.24 0.9 ]); showMap(ni,CP,1,0.5,0.1)
axes('Position',[1-1/4*(3)+0.0 0.01 0.24 0.9 ]); showMap(ni,CQ,1,4,0.5);
axes('Position',[1-1/4*(2)+0.0 0.01 0.24 0.9 ]); showMap(ni,CST,1,2,0.25);
axes('Position',[1-1/4*(1)+0.0 0.01 0.24 0.9 ]); showMap(ni,CT,1,2,0.25);

if ni >0
	saveas(gcf,[dir 'SurfAmplpots ' num2str(ni) ' ' fn]); 
else
	saveas(gcf,[dir 'SurfAmplpots_Surf ' fn]); 	
end
end
%%
if doQRSArea
hf=figure(2032);set(hf,'PaperPositionMode','auto')
QRSint=sum((PHI(:,av:floor(av+max(depV)))),2);
set(gcf,'position',[520   100   1011 406  ]);clf; 
colormap(AS.AVO);%(204:end,:))
showMapST(ni,QRSint,1,120,25)
if ni >0
	saveas(gcf,[dir 'QRSint ' num2str(ni) ' ' fn]); 
else
	saveas(gcf,[dir 'QRSint_Surf ' fn]); 	
end
end
%%
if doST
hf=figure(2033);set(hf,'PaperPositionMode','auto')
set(gcf,'position',[520   100   1011 406  ]);clf; 
colormap(AS.AVO)
showMapST(ni,CST,1,2,0.25)
if ni >0
	saveas(gcf,[dir 'SurfST ' num2str(ni) ' ' fn]); 
else
	saveas(gcf,[dir 'SurfST_Surf ' fn]); 	
end
end


