% function invqtri(qtrimode, figname, GEOM, data, scanmode)
%
% Two kind of figures can be created and saved;
%   1. [qtrimode = 1]: Correlation and RD values shown on heart surface
%   2. [qtrimode = 2]: Depolarization times for initial and final solution
%
% INPUT:
% qtrimode      = [1] Correlation and RD values shown on heart surface, 
%                 [2] Depolarization times for initial and final solution
% figname       = Name for figures to be saved as, inlcuding pathname. Without .png
% GEOM          = Structure created with InvInit.m.
% data          = Depends on qtrimode: 
%                 [qtrimode = 1] Contains .corsinit, .rdsinit and .measinit
%                 [qtrimode = 2] Contains .meas and .measinit
% scanmode      = mode for getSmode [default == 1]
%
% Version 1: 01-apr-2015

function invqtri(qtrimode, figname, GEOM, data, scanmode)

if ~exist('scanmode','var'), scanmode = 1; end

if qtrimode == 1,
    
    corsinit    = data.corsinit;
    rdsinit     = data.rdsinit;
    measinit    = data.measinit;
    
    for k = 1:size(corsinit,2)
        
        % visualize correlation plot:
        qtriplot;
        qtriplot('reset');          pause(0.2);
        qtriplot('horizontal 2');
        qtriplot('bgdcolor white'); pause(0.1);
        qtriplot('size 1000 800');  pause(0.1);
        
        % panel 1,1
        qtriplot(GEOM.VER,GEOM.ITRI,'surf');
        qtriplot(corsinit(:,k));
        qtriplot(GEOM.VER(measinit.foci(k),:),[],'fociinit');
        qtriplot('color green');
        
        % panel 1,2
        qtriplot(GEOM.VER,GEOM.ITRI,'surf_2');
        qtriplot(rdsinit(:,k));
        qtriplot(GEOM.VER(measinit.foci(k),:),[],'fociinit_2');
        qtriplot('color green');
        qtriplot('panel surf_2=2 1');
        qtriplot('panel fociinit_2=2 1');
        
        qtriplot('funcolor uniheat');
        qtriplot(['funscale 0 ' num2str(max(corsinit(:,k)))]);
        qtriplot('step 0.02');
        qtriplot('scale 60');
        
        pause
        
        % save figure:
        qtriplot(['png ' figname '_focus_' num2str(measinit.foci(k)) '.png 1000 1000']); pause(0.2);
        
        qtriplot('exit');
        
    end
    
elseif qtrimode == 2,
    
        meas        = data.meas;
        measinit    = data.measinit;
        
        if scanmode == 4,
            meas.depfinal   = meas.repfinal;
            measinit.dep    = measinit.rep;
        end
    
    [mindep,imindep]    = min(meas.depfinal);                                       % first site of depolarization
    idepearly           = find(meas.depfinal<=(mindep+2) & meas.depfinal>(mindep)); % sites within 2 ms of first site
    
    % Determine difference between earliest endo and epicardiac depolarization times
    indendo             = find(meas.depfinal == min(meas.depfinal(GEOM.endoVER==1)));
    int_endo            = find(GEOM.endoVER(indendo) == 1);
    measendo            = indendo(int_endo);
    
    indepi              = find(meas.depfinal == min(meas.depfinal(GEOM.endoVER~=1)));
    int_epi             = find(GEOM.endoVER(indepi) == 0);
    measepi             = indepi(int_epi);
    
    diffendoepi         = meas.depfinal(measendo) - meas.depfinal(measepi);
    
    indendo_init        = find(measinit.dep == min(measinit.dep(GEOM.endoVER==1)));
    int_endo            = find(GEOM.endoVER(indendo_init) == 1);
    measendo_init       = indendo_init(int_endo);
    
    indepi_init        = find(measinit.dep == min(measinit.dep(GEOM.endoVER~=1)));
    int_epi            = find(GEOM.endoVER(indepi_init) == 0);
    measepi_init       = indepi_init(int_epi);
    
    diffendoepi_init    = measinit.dep(measendo_init) - measinit.dep(measepi_init);
    
    %% visualize initial estimation and final solution inverse with qtriplot [DEPOLARIZATION]:
    aptg    = 'angle 160 20 10';
    sz      = 'step 5';
    
    qtriplot('reset');
    qtriplot('horizontal 2');
    qtriplot('vertical 2');
    qtriplot('bgdcolor white'); pause(0.1);
    qtriplot('size 1000 1000'); pause(0.1);
    
    % panel 1,1
    qtriplot(GEOM.VER,GEOM.ITRI,'init_1');
    qtriplot(measinit.dep);
    qtriplot(aptg)
    qtriplot(GEOM.VER(measendo_init,:),[],'fociendoinit_1');
    qtriplot('color cyan');
    qtriplot(aptg)
    qtriplot(GEOM.VER(measepi_init,:),[],'fociepiinit_1');
    qtriplot('color cyan');
    qtriplot(aptg)
    qtriplot(GEOM.VER(measinit.foci,:),[],'fociinit_1');
    qtriplot('color white');
    qtriplot(aptg)
    
    % panel 1,2
    qtriplot(GEOM.VER,GEOM.ITRI,'init_2');
    qtriplot(measinit.dep);
    qtriplot(GEOM.VER(measendo_init,:),[],'fociendoinit_2');
    qtriplot('color cyan');
    qtriplot(GEOM.VER(measepi_init,:),[],'fociepiinit_2');
    qtriplot('color cyan');
    qtriplot(GEOM.VER(measinit.foci,:),[],'fociinit_2');
    qtriplot('color white');
    qtriplot('panel init_2=1 2');
    qtriplot('panel fociendoinit_2=1 2');
    qtriplot('panel fociepiinit_2=1 2');
    qtriplot('panel fociinit_2=1 2');
    
    % panel 2,1
    qtriplot(GEOM.VER,GEOM.ITRI,'final_1');
    qtriplot(meas.depfinal);
    qtriplot(aptg)
    qtriplot(GEOM.VER(measendo,:),[],'fociendo_1');
    qtriplot('color cyan');
    qtriplot(aptg)
    qtriplot(GEOM.VER(measepi,:),[],'fociepi_1');
    qtriplot('color cyan');
    qtriplot(aptg)
    qtriplot(GEOM.VER(imindep,:),[],'focusfinal_1');
    qtriplot('color green');
    qtriplot(aptg)
    qtriplot(GEOM.VER(idepearly,:),[],'focifinal_1');
    qtriplot('color white');
    qtriplot(aptg)
    qtriplot('panel final_1=2 1');
    qtriplot('panel focifinal_1=2 1');
    qtriplot('panel fociendo_1=2 1');
    qtriplot('panel fociepi_1=2 1');
    qtriplot('panel focusfinal_1=2 1');
    
    
    % panel 2,2
    qtriplot(GEOM.VER,GEOM.ITRI,'final');
    qtriplot(meas.depfinal);
    qtriplot(GEOM.VER(measendo,:),[],'fociendo');
    qtriplot('color cyan');
    qtriplot(GEOM.VER(measepi,:),[],'fociepi');
    qtriplot('color cyan');
    qtriplot(GEOM.VER(imindep,:),[],'focusfinal');
    qtriplot('color green');
    qtriplot(GEOM.VER(idepearly,:),[],'focifinal');
    qtriplot('color white');
    qtriplot('panel final=2 2');
    qtriplot('panel focifinal=2 2');
    qtriplot('panel fociendo=2 2');
    qtriplot('panel fociepi=2 2');
    qtriplot('panel focusfinal=2 2');
    
    % color scale
    qtriplot('funcolor tim');
    qtriplot('funscale autocol');
    qtriplot(sz);
    
    % show values:
    qtriplot(['text 0.10 0.98 RD    = ' num2str(measinit.rd) '~0.5']); pause(0.1);
    qtriplot(['text 0.10 0.96 COR   = ' num2str(measinit.cor) '~0.5']); pause(0.1);
    qtriplot(['text 0.10 0.94 FOCUS = ' num2str(measinit.foci) '~0.5']); pause(0.1);
    qtriplot(['text 0.25 0.98 ENDO-EPI = ' num2str(diffendoepi_init) ' ms ~0.5']); pause(0.1);
    
    qtriplot(['text 0.50 0.98 RD   = ' num2str(meas.rdfinal) '~0.5']); pause(0.1);
    qtriplot(['text 0.50 0.96 COR  = ' num2str(meas.corfinal) '~0.5']); pause(0.1);
    qtriplot(['text 0.50 0.94 ITER = ' num2str(meas.iterfinal) '~0.5']); pause(0.1);
    qtriplot(['text 0.65 0.98 FOCUS = ' num2str(imindep) '~0.5']); pause(0.1);
    
    if size(idepearly,1) > 5, idepearly = idepearly(1:5); end
    qtriplot(['text 0.65 0.96 < 2ms = ' num2str(idepearly') '~0.5']); pause(0.1);
    qtriplot(['text 0.65 0.94 ENDO-EPI = ' num2str(diffendoepi) ' ms ~0.5']); pause(0.1);
    
    qtriplot('scale 70');
    
    pause
    
    % save figure:
    qtriplot(['png ' figname '.png 1000 1000']); pause(0.2);
    
    qtriplot('exit');
    
end