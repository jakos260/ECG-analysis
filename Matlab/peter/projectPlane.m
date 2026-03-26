function projVer = projectPlane(planeVectorIn, planePosIn, VER)

planeVector = planeVectorIn ./ norm(planeVectorIn);
planePos = - dot( planePosIn, planeVector );

projVer  = VER - ( dot(VER', (ones(size(VER,1),1) * planeVector)' )' + planePos) * planeVector;