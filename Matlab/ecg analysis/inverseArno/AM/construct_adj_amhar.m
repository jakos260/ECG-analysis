
% load ventricles:
[heart.VER, heart.ITRI] = loadtri('amhar.tri');

%% construct adjecency matrices:
tic
[heart.ADJsurf, heart.DISTsurf, heart.PATH] = graphdist(heart.ITRI,heart.VER,4);
[heart.ADJ, heart.DIST, heart.PATH3D]       = graphdist(heart.ITRI,heart.VER,3);
toc
%%
% load endocardial surface geometries:
[heart.LVER,heart.LITRI] = loadtri('amlho.tri'); 
[heart.RVER,heart.RITRI] = loadtri('amrho.tri');

%  construct anisotropic matrices:
heart.anisotropyRatio    = 2;                                                   % is anisotropy for conduction velocity not conductivity!!
heart.ADJ2W              = calcAnisADJ(heart,heart.anisotropyRatio);            % include anisotropy in adjecency matrix
heart.DIST2W             = graphdist(heart.ADJ2W);                              % include anisotropy in distance matrix

%%

save('adj_amhar.mat','heart')