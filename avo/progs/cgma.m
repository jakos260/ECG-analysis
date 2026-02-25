% cgma.m
% computes the transfer between UDL nodes on a triangulated surface (VER,ITRI) 
% and observation points OBS;
% quick and dirty version: assumes observation points to be remote
% A. van Oosterom; 20041004
function GMA=cgma(VER,ITRI,OBS)
nobs=size(OBS,1);
nver=size(VER,1);
ntri=size(ITRI,1);
GMA=zeros(nobs,nver);
for i=1:nobs;
    solids=solida(VER,ITRI,OBS(i,:));
    GMA(i,ITRI(:,1))=GMA(i,ITRI(:,1))+solids;
    GMA(i,ITRI(:,2))=GMA(i,ITRI(:,2))+solids;
    GMA(i,ITRI(:,3))=GMA(i,ITRI(:,3))+solids;
end

GMA=GMA/(6*pi);