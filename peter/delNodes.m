function [VER,ITRI]=delNodes(VER,ITRI)

keep = unique(ITRI);
if ( length(keep) ~= length(VER) )
    remVer= 1:length(VER);
    keepVer= 1:length(VER);
    remVer(keep)=[];
    keepVer(remVer)=[];
    VER(remVer,:)=[];
    I = ITRI;
    for i=1:size(VER,1)
        I(I==keepVer(i)) = i;
    end
    ITRI = I;
end







