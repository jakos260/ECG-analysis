

%% visualize initial estimation and final solution inverse with qtriplot [DEPOLARIZATION]:
aptg = 'angle 160 20 10';

qtriplot('reset');
qtriplot('horizontal 2');
qtriplot('vertical 2');
qtriplot('bgdcolor white'); pause(0.1);
qtriplot('size 1000 1000'); pause(0.1);

% panel 1,1
qtriplot(GEOM.VER,GEOM.ITRI,'init_1');
qtriplot(measinit.dep);
qtriplot(aptg)
qtriplot(GEOM.VER(measinit.foci,:),[],'fociinit_1');
qtriplot('color white');
qtriplot(aptg)

% panel 1,2
qtriplot(GEOM.VER,GEOM.ITRI,'init_2');
qtriplot(measinit.dep);
qtriplot(GEOM.VER(measinit.foci,:),[],'fociinit_2');
qtriplot('color white');
qtriplot('panel init_2=1 2');
qtriplot('panel fociinit_2=1 2');

% panel 2,1
qtriplot(GEOM.VER,GEOM.ITRI,'final_1');
qtriplot(meas.depfinal);
qtriplot(aptg)
qtriplot(GEOM.VER(imindep,:),[],'focusfinal_1');
qtriplot('color green');
qtriplot(aptg)
qtriplot(GEOM.VER(idepearly,:),[],'focifinal_1');
qtriplot('color white');
qtriplot(aptg)
qtriplot('panel final_1=2 1');
qtriplot('panel focusfinal_1=2 1');
qtriplot('panel focifinal_1=2 1');

% panel 2,2
qtriplot(GEOM.VER,GEOM.ITRI,'final');
qtriplot(meas.depfinal);
qtriplot(GEOM.VER(imindep,:),[],'focusfinal');
qtriplot('color green');
qtriplot(GEOM.VER(idepearly,:),[],'focifinal');
qtriplot('color white');
qtriplot('panel final=2 2');
qtriplot('panel focusfinal=2 2');
qtriplot('panel focifinal=2 2');

% color scale
qtriplot('funcolor tim');
qtriplot('funscale autocol');

% show values:
qtriplot(['text 0.1 0.98 RD    = ' num2str(measinit.rd) '~0.5']); pause(0.1);
qtriplot(['text 0.1 0.96 COR   = ' num2str(measinit.cor) '~0.5']); pause(0.1);
qtriplot(['text 0.1 0.94 FOCUS = ' num2str(measinit.foci) '~0.5']); pause(0.1);

qtriplot(['text 0.5 0.98 RD   = ' num2str(meas.rdfinal) '~0.5']); pause(0.1);
qtriplot(['text 0.5 0.96 COR  = ' num2str(meas.corfinal) '~0.5']); pause(0.1);
qtriplot(['text 0.5 0.94 ITER = ' num2str(meas.iterfinal) '~0.5']); pause(0.1);
qtriplot(['text 0.7 0.98 FOCUS = ' num2str(imindep) '~0.5']); pause(0.1);
qtriplot(['text 0.7 0.96 EARLY = ' num2str(idepearly') '~0.5']); pause(0.1);

qtriplot(['text 0.1 0.06 mudep = ' num2str(mudep) '~0.3']); pause(0.1);
qtriplot(['text 0.1 0.05 murep = ' num2str(murep) '~0.3']); pause(0.1);
qtriplot(['text 0.1 0.04 INIT ANIS = ' num2str(initanis) '~0.3']); pause(0.1);
qtriplot(['text 0.1 0.03 ANIS = ' num2str(anis) '~0.3']); pause(0.1);
qtriplot(['text 0.1 0.02 INIT VEL = ' num2str(initialvelocity) '~0.3']); pause(0.1);

pause

% save figure:
qtriplot(['png ' S.dirout '/figures/focus_' fnameout_fig  '.png 1000 1000']); pause(0.1);

qtriplot('exit');