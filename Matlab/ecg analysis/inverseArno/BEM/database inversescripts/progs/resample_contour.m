% resample_contour.m
% 
% NEWC=resample_contour(CONT,ncrec,fig);
% place ncrec nodes on a 3D contour; if the contour is closed, CONT(1,:)
% should be equal to CONT(size(CONT,1),:)
% linear interpolation between subsequent nodes of the contour
% A. van Oosterom; 20110904
% see also: NEWC=fourier_contours(CONT,kmax,fig,ncrec)

function NEWC=resample_contour(CONT,ncrec,fig)

nct=size(CONT,1);
closed=0;

if sum(CONT(1,:)==CONT(nct,:),2)==3,
    closed=1;
end

if nargin>2,
      figure(fig)
      hold on
      plot3(CONT(:,1),CONT(:,2),CONT(:,3))
      plot3(CONT(:,1),CONT(:,2),CONT(:,3),'+')
      grid on 
end

% s = (rectified) path length along the contour
s=[0; cumsum(norm3d(CONT(2:nct,:)-CONT(1:nct-1,:)))];
smax=s(nct);

% resample
NEWC=zeros(ncrec,3);
NEWC(1,:)=CONT(1,:);
NEWC(ncrec,:)=CONT(nct,:);

if nargin>2,
    plot3(NEWC(1,1),NEWC(1,2),NEWC(1,3),'k*','linewidth',1.5)
end
dels=smax/(ncrec-1);
for j=2:ncrec-1,
    t=(j-1)*dels;
    k=find(s<=t);
    k=max(k);
    %[s(k) t s(k+1) ]
    f=(t-s(k))/(s(k+1)-s(k));
    NEWC(j,:)=(1-f)*CONT(k,:)+f*CONT(k+1,:);
    if nargin>2,
       plot3(NEWC(j,1),NEWC(j,2),NEWC(j,3),'r*','linewidth',1.5)
    end
     
end

if nargin>2,
    plot3(NEWC(ncrec,1),NEWC(ncrec,2),NEWC(ncrec,3),'r*','linewidth',1.5)
end
end






