function [VER,ITRI] = cleanTRIs(VER,ITRI)

for i=1:length(VER)
    a=find(VER(i+1:end,1)==VER(i,1) & VER(i+1:end,2)==VER(i,2) & VER(i+1:end,3)==VER(i,3));
    if ~isempty(a)
        for j=1:length(a)
            ITRI(ITRI==a(j)+i)=i;
        end
    end
end
dels = zeros(length(ITRI),1);
for i=1:length(ITRI)
    if ITRI(i,1)==ITRI(i,2) || ITRI(i,1)==ITRI(i,3) || ITRI(i,2)==ITRI(i,3)
        dels(i)=1;
    end
end
ITRI(dels==1,:)=[];

i=1;
while i <= length(VER)
    if isempty(find(ITRI==i,1))
        [VER,ITRI]=delNode(VER,ITRI,i);
    else
        i=i+1;
    end
end
savetri('test.tri',VER,ITRI);