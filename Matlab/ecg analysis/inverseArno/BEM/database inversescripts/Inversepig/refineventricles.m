if 0
    clearvars
    [VER,ITRI]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_ventricles.tri');
    [VERr,ITRIr]=refine4(VER,ITRI);
    typ=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_ventricles.typ');
    
    [VERcr,ITRIcr]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_rcav.tri');
    [VERcrr,ITRIcrr]=refine4(VERcr,ITRIcr);
    savetri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09r/Pig09r_rcav.tri',VERcrr,ITRIcrr);
    [VERcl,ITRIcl]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_lcav.tri');
    [VERclr,ITRIclr]=refine4(VERcl,ITRIcl);
    savetri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09r/Pig09r_lcav.tri',VERclr,ITRIclr);
    
    [VERa,ITRIa]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_atria.tri');
    [VERar,ITRIar]=refine4(VERa,ITRIa);
    savetri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09r/Pig09r_atria.tri',VERar,ITRIar);
end
nVER1=size(VER,1);
nVER2=size(VERr,1);


typ2=NaN(nVER2,1);
typ2(1:nVER1,1)=typ;
for i=(nVER1+1):nVER2
    iver=find(ITRIr(:,1)==i | ITRIr(:,2)==i | ITRIr(:,3)==i);
    in=ITRIr(iver,:);
    % find edge with both <
    in=in(in<=nVER1);
    if all(typ(in)>3) && all(typ(in)==typ(in(1)))
    typ2(i)=typ(in(1));
    else
        typ2(i)=min(typ(in));
    end
    
end

if 0
    saveasciint('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09r/Pig09r_ventricles.typ',typ2);
    savetri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09r/Pig09r_ventricles.tri',VERr,ITRIr);
end

if 0 % test
    typtest=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09r/Pig09r_ventricles.typ');
    [VERt,ITRIt]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09r/Pig09r_ventricles.tri');
    qtriplot({{'delete *'}
        {VERt,ITRIt,'test'}
        {typtest,'test'}});
end

if 0 % test saving
   savetri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/test/test_ventricles.tri', VER, ITRI);
   savetri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/test/test_atria.tri', VERa, ITRIa);
   savetri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/test/test_lcav.tri', VERcl, ITRIcl);
   savetri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/test/test_rcav.tri', VERcr, ITRIcr);
    
end
        
% [VERtest,ITRItest]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09r/Pig09r_atria.tri');