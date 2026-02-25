% show_nimelecs(VER,mode,mcolor)
% show electrode position of the nim64 lead system; 
% mode=-9 show standard 12; mode =9 includes labeling
% mode=-64 shows all 64 lead electrodes; mode=64 includes the labels: 1:64
% default: mode=-64
% mcolor=marker color
% a. van oosterom; 20080708

function show_nimelecs(VER,mode,mcolor)
if nargin==1, mode=-64; end
if abs(mode)~=9 & abs(mode)~=64, mode=-64; end
if nargin<3, mcolor='w'; end
hold on

if abs(mode)==9,
   % llabels !!!!!! tuned to nim64_leads
   L9 = [
   VER(1,:)  % VR 
   VER(2,:)  % VL
   VER(87,:) % VF
   VER(19,:) % V1
   VER(26,:) % V2
   (VER(33,:)+VER(34,:))/1.95; % V3  
%  NBNB: there was no electrode at this location; the local potential
%  as based on measured data could be assigned to be:   PHI((VER(33,:))+PHI(VER(34,:))/2,
%  or mean(PHI(VER[ 26 27 33 34 40 41],:),:)
%
   VER(41,:) % V4
   VER(48,:) % V5
   VER(54,:) % V6 
   ]*1.02;
   
   mark9=line(L9(:,1),L9(:,2),L9(:,3),'LineStyle','none','Marker','o',...
   'MarkerSize',6,'MarkerFaceColor',mcolor,'MarkerEdgeColor',mcolor);
   name = {'VR','VL','VF','V1','V2','V3','V4','V5','V6'};
   fact = 1.15;
   
   if mode>0,
      for i=1:9
	      mark9labs= text(L9(i,1)*fact,L9(i,2)*fact,L9(i,3)*fact+0.03,name{i});
	      set(mark9labs,'FontSize',10);
      end
  end
end

if abs(mode)==64,
    mark64=line(VER(1:64,1),VER(1:64,2),VER(1:64,3),'LineStyle','none','Marker','o',...
   'MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','w');
   if mode>0,
      fact = 1.02;
      for i=1:64
	      mark9labs= text(VER(i,1)*fact,VER(i,2)*fact,VER(i,3)*fact+0.005,num2str(i));
	      set(mark9labs,'FontSize',10);
      end
  end
end