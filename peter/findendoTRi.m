
function ind=findendoTRi(VER,ITRI,EVER)

mask=zeros(length(ITRI),1);
ind=[1:length(ITRI)];

for i=1:length(ITRI)
	epivers=VER(ITRI(i,:),:);
	if ~isempty(find(EVER(:,1)==epivers(1,1) & EVER(:,2)==epivers(1,2) & EVER(:,3)==epivers(1,3)) &...
					find(EVER(:,1)==epivers(2,1) & EVER(:,2)==epivers(2,2) & EVER(:,3)==epivers(2,3)) &...  
					find(EVER(:,1)==epivers(3,1) & EVER(:,2)==epivers(3,2) & EVER(:,3)==epivers(3,3)))
		mask(i)=1;
	end
end
ind=ind(mask==1);