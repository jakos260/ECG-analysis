function ind=findVERi(VER1,VER2)

ind=[];
for i=1:length(VER1)
    if ~isempty(find(VER1(i,1)==VER2(:,1) & VER1(i,2)==VER2(:,2) & VER1(i,3)==VER2(:,3),1))
        ind=[ind; i];
    end
end