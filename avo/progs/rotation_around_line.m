% rotation_around_line.m
% ROTS=rotation_around_line(PNTS,alpha,begpoint,endpoint)
% PNTS: matrix of n points in 3D: size(n,3)
% alpha: desired rotation (radians)
% begpoint and endpoint define the axis around which the rotation is
% performed
% A. van Oosterom; 20140819

    function ROTS=rotation_around_line(PNTS,alpha,begpoint,endpoint)
    np=size(PNTS,1);
    as=endpoint-begpoint;
    nas=as/norm(as);
    
    NAS=ones(np,1)*nas;
    NAXIS1=ones(np,1)*begpoint;
    % projections of PNTS on the line (axis)
    PROs=NAXIS1 + diag(sum((PNTS-NAXIS1).*NAS,2))* NAS;
    N=cross(NAS,(PNTS-PROs));
    ROTS=PROs + cos(alpha)*(PNTS-PROs) + sin(alpha)*N;   
   
        
        
   
    

