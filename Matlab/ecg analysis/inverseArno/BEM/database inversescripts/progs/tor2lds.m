function T=tor2lds(leadspecs,surface)
% computes the transfer from a scalar function defined at 
% nn nodes of a triangulated surface 
% to nl positions specified by lambda and mu value of the 
% triangle that carries the lead
% 2003-04-18

% METHOD:
% biliniar interpolation
LEADS=loadmat(leadspecs);
[nl,jdum]=size(LEADS);
[VER,ITRI]=loadtri(surface);
[nn,jdum]=size(VER);
T=zeros(nl,nn);
for i=1:nl,
   ind=ITRI(LEADS(i,1),:);
   T(1,ind(1))=1-LEADS(i,2)-LEADS(i,3);
   T(i,ind(2))=LEADS(i,2);
   T(i,ind(3))=LEADS(i,3);
   end   
