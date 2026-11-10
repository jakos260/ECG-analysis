function [AMA_A,AMA_V] = getAMA(DATA,ver)
    if ~isfield(DATA, 'GEOM') || ~isfield(DATA.GEOM, 'thorax') || ~isfield(DATA.GEOM.thorax, 'VER')
        AMA_A = []; AMA_V = []; return;
    end
    
    meanTh = mean(DATA.GEOM.thorax.VER);
    if isfield(DATA,'ATRIA') && isfield(DATA.ATRIA,'THORAX')
        AMA_A = zeros(size(ver,1),size(DATA.ATRIA.THORAX,2));
    else
        AMA_A = [];
    end
    
    if isfield(DATA, 'VENTR') && isfield(DATA.VENTR, 'THORAX')
        AMA_V = zeros(size(ver,1),size(DATA.VENTR.THORAX,2));
        for i=1:size(ver,1)
            TRIS= linetris(DATA.GEOM.thorax.VER,DATA.GEOM.thorax.ITRI,ver(i,:),[meanTh(1) meanTh(2) ver(i,3)]);
            TRIS(abs(TRIS(:,end))> min(abs(TRIS(:,end))),:)=[];
            itri = DATA.GEOM.thorax.ITRI(TRIS(1),:);
            if isfield(DATA,'ATRIA') && isfield(DATA.ATRIA,'THORAX')
                AMA_A(i,:) = DATA.ATRIA.THORAX(itri(1),:) + ...
                    (DATA.ATRIA.THORAX(itri(2),:) - DATA.ATRIA.THORAX(itri(1),:)) * TRIS(1,3) + ...
                    (DATA.ATRIA.THORAX(itri(3),:) - DATA.ATRIA.THORAX(itri(1),:)) * TRIS(1,4);
            end
            AMA_V(i,:) = DATA.VENTR.THORAX(itri(1),:) + ...
                (DATA.VENTR.THORAX(itri(2),:) - DATA.VENTR.THORAX(itri(1),:)) * TRIS(1,3) + ...
                (DATA.VENTR.THORAX(itri(3),:) - DATA.VENTR.THORAX(itri(1),:)) * TRIS(1,4);
        end
    else
        AMA_V = [];
    end
end