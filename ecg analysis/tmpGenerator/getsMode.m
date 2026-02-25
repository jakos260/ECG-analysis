classdef getsMode
    enumeration
        Exponent            % -> S=1./(1+exp(-p(1)*TDEP))
        LineraST            % -> S=upslope+downslope with linear ST_section control
        CosST               % -> S=upslope+downslope with cosinus ST_section control
        LinST_rept_dept_ang % -> S=upslope+downslope+linear ST_section | p(rept, dept, st_ang)
        LinST_rept_ang_dept % -> S=upslope+downslope+linear ST_section | p(rept, st_ang, dept)
        NoReptCorrection    % -> S without rept correction | p(rept, st_ang, dept)
        Exp_Spline          % -> S=uplope + spline
        UpslopeDownslope    % -> S=upslope+downslope
    end
end