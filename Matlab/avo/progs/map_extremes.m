% map_extremes
% display sites of precomputed local extremes : loc max and loc_min

if ~exist('locmax'),
    'locmax and or locsaddle  and or locmin should have been defined'
    return
end


if exist('locmax'),
    plot3(VER(locmax,1),VER(locmax,2),VER(locmax,3),'r*','linewidth',1.5)
end

if exist('locsaddle'),
    plot3(VER(locsaddle,1),VER(locsaddle,2),VER(locsaddle,3),'w*','linewidth',1.5)
end

if exist('locmin')
    plot3(VER(locmin,1),VER(locmin,2),VER(locmin,3),'b*','linewidth',1.5)
    
end
