function PHI= phifrommeas(GEOM,meas,lpass)

t=0:min(GEOM.SPECS.endtwave,800) - GEOM.SPECS.onsetqrs;
T=ones(length(GEOM.VER),1)*t;

PHI = lowpassma(GEOM.AMA * getSmode(T,meas.depfinal,meas.repfinal,GEOM.SPECS,4),lpass);