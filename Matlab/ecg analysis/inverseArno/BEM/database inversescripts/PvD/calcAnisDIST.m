
function GEOM=calcAnisDIST(GEOM,anisotropyRatio)

if anisotropyRatio==1
	GEOM.DIST2W=GEOM.DIST;
	GEOM.ADJ2W=GEOM.ADJ;
	return;
end
GEOM.ADJ2W=calcAnisADJ(GEOM,anisotropyRatio);
tic;GEOM.DIST2W=graphdist(GEOM.ADJ2W);toc
