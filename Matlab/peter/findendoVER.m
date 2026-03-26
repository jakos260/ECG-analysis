
function ind=findendoVER(VER,ITRI,EVER)

mask=zeros(length(VER),1);
ind=[1:length(VER)];

for i=1:length(EVER)
	if ~isempty(find(EVER(:,1)==VER(i,1) & EVER(:,2)==VER(i,2) & EVER(:,3)==VER(i,3)))
		mask(i)=1;
	end
end
ind=ind(mask==1);