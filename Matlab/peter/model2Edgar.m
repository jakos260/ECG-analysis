function model2Edgar(dirIn,dirOut)

if ~exist(dirout)
    mkdir(dirout)
end

saveheart2Edgar([dirIn '_ventricles.tri'],[dirOut '_ventricles.tri'])
savetri2Edgar([dirIn '_lad.tri'],[dirOut '_lad.tri'])
savetri2Edgar([dirIn '_rca.tri'],[dirOut '_rca.tri'])
savetri2Edgar([dirIn '_lcx.tri'],[dirOut '_lcx.tri'])

saveheart2Edgar([dirIn '_atria.tri'],[dirOut '_atria.tri'])

savetri2Edgar([dirIn '_lcav.tri'],[dirOut '_lcav.tri'])
savetri2Edgar([dirIn '_rcav.tri'],[dirOut '_rcav.tri'])
savetri2Edgar([dirIn '_llung.tri'],[dirOut '_llung.tri'])
savetri2Edgar([dirIn '_rlung.tri'],[dirOut '_rlung.tri'])

savetri2Edgar([dirIn '_thorax.tri'],[dirOut '_thorax.tri'])





function saveheart2Edgar(fnIn,fnout)

[VER,ITRI]= loadtri(fnIn);
geom.pts = VER;
geom.fac = ITRI;

save([fnout '.mat'],'geom');


function savetri2Edgar(fnIn,fnout)

[VER,ITRI]= loadtri(fnIn);
geom.pts = VER;
geom.fac = ITRI;
save([fnout '.mat'],'geom');