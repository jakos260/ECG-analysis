% lamu2xyz.m
function XYZ=lamu2xyz(LAMU,ITRI,VER)

% for a point specified by lambda and mu on a triangle:
% find its xyz co-ordinates

dim=size(LAMU);
npnts=dim(1);
dim=size(ITRI);
ntri=dim(1);
dim=size(VER);
nver=dim(1);
% for each point of the npnts points specified in 
% LAMU by lambda and mu values: find its xyz values
for n=1:npnts,
itr=LAMU(n);
i=zeros(1,3);
i=ITRI(itr,:);
XYZ(n,:)=VER(i(1),:)+(VER(i(2),:)-VER(i(1),:))*LAMU(n,2)+(VER(i(3),:)-VER(i(1),:))*LAMU(n,3);
end

