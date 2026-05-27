m       = mean(GEOM.VER);                               % determine center of ventricles mesh
T       = linetris(GEOM.TVER,GEOM.TITRI,m,m+[1,0,0]);   % find intersection mesh with line through center in x-direction
i1      = find(T(:,2)<0);i2=find(T(:,2)>0);
v1      = m+T(i1,5)*[1 0 0];                            % calculate intersection line with triangle v1
v2      = m+T(i2,5)*[1 0 0];                            % calculate intersection line with triangle v2


leads   = []; % --> deze word steeds verder uitgebreid!!!
m       = v1 + 0.4*(v2-v1);

dV      = [GEOM.TVER(:,1)-m(1) GEOM.TVER(:,2)-m(2) GEOM.TVER(:,3)-m(3)];
nv      = norm3d(dV);

T       = linetris(GEOM.TVER,GEOM.TITRI,m,m+[0,1,0]);
i       = GEOM.TITRI(T(1,1),:); if T(1,3)>.5, v = i(2); elseif T(1,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];
i       = GEOM.TITRI(T(2,1),:); if T(2,3)>.5, v = i(2); elseif T(2,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];

T       = linetris(GEOM.TVER,GEOM.TITRI,m,m+[1,0,0]);
i       = GEOM.TITRI(T(1,1),:); if T(1,3)>.5, v = i(2); elseif T(1,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];
i       = GEOM.TITRI(T(2,1),:); if T(2,3)>.5, v = i(2); elseif T(2,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];

T       = linetris(GEOM.TVER,GEOM.TITRI,m,m+[1,1,0]);
i       = GEOM.TITRI(T(1,1),:); if T(1,3)>.5, v = i(2); elseif T(1,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];
i       = GEOM.TITRI(T(2,1),:); if T(2,3)>.5, v = i(2); elseif T(2,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];

T       = linetris(GEOM.TVER,GEOM.TITRI,m,m+[1,-1,0]);
i       = GEOM.TITRI(T(1,1),:); if T(1,3)>.5, v = i(2); elseif T(1,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];
i       = GEOM.TITRI(T(2,1),:); if T(2,3)>.5, v = i(2); elseif T(2,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];

% ==============================
m       = v1 + 0.25*(v2-v1);

T       = linetris(GEOM.TVER,GEOM.TITRI,m,m+[0,0,1]);
i       = GEOM.TITRI(T(1,1),:); if T(1,3)>.5, v = i(2); elseif T(1,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];
i       = GEOM.TITRI(T(2,1),:); if T(2,3)>.5, v = i(2); elseif T(2,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];

T       = linetris(GEOM.TVER,GEOM.TITRI,m,m+[0,-1,1]);
i       = GEOM.TITRI(T(1,1),:); if T(1,3)>.5, v = i(2); elseif T(1,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];
i       = GEOM.TITRI(T(2,1),:); if T(2,3)>.5, v = i(2); elseif T(2,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];

T       = linetris(GEOM.TVER,GEOM.TITRI,m,m+[0,1,1]);
i       = GEOM.TITRI(T(1,1),:); if T(1,3)>.5, v = i(2); elseif T(1,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];
i       = GEOM.TITRI(T(2,1),:); if T(2,3)>.5, v = i(2); elseif T(2,4)>.5, v = i(3); else v = i(1); end; leads   = [leads; v];


%
leads   = [];
m040 = v1 + 0.4*(v2-v1);
m025 = v1 + 0.25*(v2-v1);

dV      = [GEOM.TVER(:,1)-m040(1) GEOM.TVER(:,2)-m040(2) GEOM.TVER(:,3)-m040(3)];
nv      = norm3d(dV);

mmat = [0,1,0; 1,0,0; 1,1,0; 1,-1,0; 0,0,1; 0,-1,1; 0,1,1];

for mvar = 1:7,
    if mvar < 5, m = m040; else m = m025; end
    T = linetris(GEOM.TVER,GEOM.TITRI,m,m+mmat(mvar,:));
    for kvar = 1:2,
        i = GEOM.TITRI(T(kvar,1),:);
        if T(kvar,3)>.5, v = i(2); elseif T(kvar,4)>.5, v = i(3); else v = i(1); end;
        leads   = [leads; v];
    end
end
