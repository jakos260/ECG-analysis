
function deltaSTMap(ni,ECG,PHI,av,depV,outdir ,fn)

global TORSO
global AS
global VENTR

if ni >0
	PHI=PHI-ones(length(TORSO.TVER),1)*PHI(ni,:);
	ECG=ECG-ones(length(TORSO.TVER),1)*ECG(ni,:);
end
inforb=GetAmplInfos(ECG,max(depV)+av );
infor1=GetAmplInfos(PHI,max(depV)+av );
CP=infor1(:,1)-inforb(:,1);
CQ=infor1(:,2)-inforb(:,2);
CT=infor1(:,3)-inforb(:,3);
CST=infor1(:,4)-inforb(:,4);


doamplmap=0;
doQRSarea=1;
dodeltaST=1;
%%


if doamplmap
hf=figure(2034);set(hf,'PaperPositionMode','auto')
set(gcf,'position',[520   100   1011 406  ]);clf; 
colormap(AS.AVO)

% axes('Position',[1-1/3*(3)+0.0 0.01 0.3 0.9 ]); showMap(ni,CP,1,0.4)
% axes('Position',[1-1/3*(2)+0.0 0.01 0.3 0.9 ]); showMap(ni,CQ,1,4)
% axes('Position',[1-1/3*(1)+0.0 0.01 0.3 0.9 ]); showMap(ni,CT,1,0.8)
axes('Position',[1-1/2*(2)+0.0 0.01 0.45 0.9 ]); showMap(ni,CQ,1,4,1)
axes('Position',[1-1/2*(1)+0.0 0.01 0.45 0.9 ]); showMap(ni,CT,1,1,0.1)
if ~isempty(fn)
	if ~isempty(strfind(fn,'.tif'))
		if ni >0
			saveas(gcf,[outdir 'DeltaAmplpots ' num2str(ni) ' ' fn]); 
		else
			saveas(gcf,[outdir 'DeltaAmplpots_Surf ' fn]); 	
		end	
	else
		if ni >0
			print('-dtiffn',[outdir 'DeltaAmplpots ' num2str(ni) ' ' fn]); 
		else
			print('-dtiffn',[outdir 'DeltaAmplpots_Surf ' fn]);
		end
	end
end
end
%%
if doQRSarea
hf=figure(2035);set(hf,'PaperPositionMode','auto')
QRSintb=sum((ECG(:,av:floor(av+max(VENTR.depV)))),2);
QRSint=sum((PHI(:,av:floor(av+max(depV)))),2);
set(gcf,'position',[520   100   1011 406  ]);clf; 
colormap(AS.AVO)
showMapST(ni,QRSint-QRSintb,1,120,10)
if ~isempty(fn)
	if isempty(strfind(fn,'.tif'))
		if ni >0
			saveas(gcf,[outdir 'deltaQRSint ' num2str(ni) ' ' fn]); 
		else
			saveas(gcf,[outdir 'deltaQRSint_Surf ' fn]); 	
		end	
	else
		if ni >0
			print('-dtiffn',[outdir 'deltaQRSint ' num2str(ni) ' ' fn]); 
		else
			print('-dtiffn',[outdir 'deltaQRSint_Surf ' fn]); 	
		end
	end
end
end
%% ======================================================
if dodeltaST
	
hf=figure(2036);set(hf,'PaperPositionMode','auto')
set(gcf,'position',[520   100   1011 406  ]);clf; 

colormap(AS.AVO)
showMapST(ni,CST,1,2,0.25)
if ~isempty(fn)
	if ~isempty(strfind(fn,'.tif'))
		if ni >0
			print('-dtiffn',[outdir 'deltaST_ ' num2str(ni) ' ' fn]); 
		else
			print('-dtiffn',[outdir 'deltaST_Surf ' fn]); 	
		end
	else
		if ni >0
			saveas(gcf,[outdir 'deltaST_ ' num2str(ni) ' ' fn]); 
		else
			saveas(gcf,[outdir 'deltaST_Surf ' fn]); 	
		end
	end
end
end