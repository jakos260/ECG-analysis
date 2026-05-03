% function [SA] = solangle(VER,ITRI)
%
% Calculation of the solid angle for triangles defined in [ITRI].
% NOTE: THIS FUNCTION HAS TO BE CHECKED!!
%
% INPUT:
% VER   = vertices of geometry in a [N x 3] matrix.
% ITRI  = triangles of geometry in a [N x 3] matrix.
%
% OUTPUT:
% SA    = solid angles for all triangles [ITRI] in a [N x 1] matrix.
%
% IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. BME-30, NO. 2, FEB 1983
% A. VAN OOSTEROM AND J. STRACKEE
%
% Arno Janssen 13-apr-2015

function [SA] = solangle(VER, ITRI)

[nrv ncolv] = size(VER);
[nri ncoli] = size(ITRI);

if ncolv ~= 3, VER  = VER';  end
if ncoli ~= 3, ITRI = ITRI'; end

nrsa    = size(ITRI,1);
SA      = zeros(nrsa,1);

for k = 1:nrsa,
    
    % select vertices that belong to triangle k:
    r1 = VER(nrsa(k,1),:);
    r2 = VER(nrsa(k,2),:);
    r3 = VER(nrsa(k,3),:);
    
    % calculate cross product for r2 & r3:
    cp  = cross(r2,r3);
    
    % calculate numerator; i.e., the volume of the parallelepiped spanned
    % by the vectors r1, r2, and r3:
    nom = dot(r1,cp);
    
    % calculate the lengths:
    n1  = norm(r1);
    n2  = norm(r2);
    n3  = norm(r3);
    
    % calculate the dot product of r1, r2 and r3:
    ip12 = dot(r1,r2);
    ip23 = dot(r2,r3);
    ip13 = dot(r1,r3);
    
    % calculate the denominator:
    den = n1*n2*n3 + ip12*n3 + ip23*n1 + ip13*n2;
    
    % obseravtion point is on triangle:
    if nom == 0 && den <= 0, SA(k) = NaN; end
    
    % determine solid angle
    SA(k) = -2.*atan2 (nom,den);
    
end