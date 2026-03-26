function convertModel2CEI(modeldir, modelname,contourfile,outdir)
% 
DATA = readGeomPeacsModel(modeldir,modelname);

Atria_geom.node = [DATA.GEOM.atria.VER(:,2) -DATA.GEOM.atria.VER(:,1) DATA.GEOM.atria.VER(:,3)];
Atria_geom.face = DATA.GEOM.atria.ITRI;
save(fullfile(outdir,'Atria'),'Atria_geom');

Ventricles_geom.node = [DATA.GEOM.ventr.VER(:,2) -DATA.GEOM.ventr.VER(:,1) DATA.GEOM.ventr.VER(:,3)];
Ventricles_geom.face = DATA.GEOM.ventr.ITRI;
save(fullfile(outdir,'Ventricles'),'Ventricles_geom');

RCav_geom.node = [DATA.GEOM.rcav.VER(:,2) -DATA.GEOM.rcav.VER(:,1) DATA.GEOM.rcav.VER(:,3)];
RCav_geom.face = DATA.GEOM.rcav.ITRI;
save(fullfile(outdir,'RCAV'),'RCav_geom');

LCav_geom.node = [DATA.GEOM.lcav.VER(:,2) -DATA.GEOM.lcav.VER(:,1) DATA.GEOM.lcav.VER(:,3)];
LCav_geom.face = DATA.GEOM.lcav.ITRI;
save(fullfile(outdir,'LCAV'),'LCav_geom');

RLung_geom.node = [DATA.GEOM.rlung.VER(:,2) -DATA.GEOM.rlung.VER(:,1) DATA.GEOM.rlung.VER(:,3)];
RLung_geom.face = DATA.GEOM.rlung.ITRI;
save(fullfile(outdir,'RLung'),'RLung_geom');

LLung_geom.node = [DATA.GEOM.llung.VER(:,2) -DATA.GEOM.llung.VER(:,1) DATA.GEOM.llung.VER(:,3)];
LLung_geom.face = DATA.GEOM.llung.ITRI;
save(fullfile(outdir,'LLung'),'LLung_geom');

Torso_geom.node = [DATA.GEOM.thorax.VER(:,2) -DATA.GEOM.thorax.VER(:,1) DATA.GEOM.thorax.VER(:,3)];
Torso_geom.face = DATA.GEOM.thorax.ITRI;
save(fullfile(outdir,'Torso'),'Torso_geom');








A = xml2struct(contourfile);
C = A(1).Contours;
names = fieldnames(C);
for i=1:length(names)   
    eval (['VER = str2num(C.' names{i} '.Points.Text);']);
    if ~isempty(VER)
        eval (['A =textscan(C.' names{i} '.DICOM_Vertices.Text,''%s %d %d'');']);
        dicoms = A{1};
        starts = A{2};
        ends = A{3};
        for j=1:length(dicoms)
            ver = VER(starts(j):ends(j),:);               
            if ~isempty(ver)
                meta.nrrd_type = 'NRRD';
                meta.dicom = dicoms{j};
                meta.tissue = names{i};
                nrrdwrite(ver,meta,fullfile(outdir,'contours',[ names{i} '_' dicoms{j} '.nrrd' ]));
            end
        end
    end
end



