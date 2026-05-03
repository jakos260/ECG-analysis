% fig_6_14
% 20060505
% used in comprehensive chapt 6
% nb: different selection of OBS compared to what is used in
% comprehensive/chapt6 scripts

clear
a=2.4; b=4;
nr=11;
nphi=65;
velocity=1  % m/s
ra=a;
rb=b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  menu
createboundandobs=1;

endo=1; % else epi

frontsandelgs=0;

sigp=1;sigs=1;

siga=6; % cavity conductivity (relative to myocard and exterior)

if createboundandobs==0,
   % compute geometry for treating boundary effects and OBS 
   ra=a;
   rb=b;
   figure(1)
   clf
   % specification of globe shaped boundary
   ntheta=61; nphibound=12;
   [VER, ITRI]=make_globe(ntheta,nphibound); % unit sphere
   
   VERb=b*VER;
   VERa=a*VER;
  
   figure(1)
   clf
   VER=VERa;
   boxsize=4.4;
   triplot
   hold on
   
   savetri('epicard.tri', VERb,ITRI);
   savetri('endocard.tri',VERa,ITRI);
   ITRIa=ITRI;
   ITRIb=ITRI;
   VER=VERa;
   % define observation points for electrograms
   %intramural
   z=(0:0.05:4.4)';
   nz=length(z);
   OBS=[zeros(nz,1) zeros(nz,1) z] ;
   
   nobs=length(OBS);
   % show field points
   plot3(OBS(1:nobs,1),OBS(1:nobs,2),OBS(1:nobs,3),'k+')
   plot3(OBS(25,1),OBS(25,2),OBS(25,3),'k*')
   savemat('ver_obs.mat',OBS)
   pause
else,
      OBS=loadmat('ver_obs.mat');
      nobs=length(OBS);
      [VERb,ITRIb]=loadtri('epicard.tri');
      [VERa,ITRIa]=loadtri('endocard.tri');
end

% identify observation points that should coincide with boundary vertices
% force these to coincide
COIN=[];
for i=1:nobs,
    for j=1:length(VERa),
        if sum(abs(VERa(j,:)-OBS(i,:)),3)<3.e-6,
           COIN=[COIN; [i j]]; break,
        end
    end
end
OBS(COIN(:,1),:)=VERa(COIN(:,2),:);

if endo==0,
      'epicardial stumulus'
      pnt1=[0 b 0];
else,
      'endocardial stimulus'
      pnt1=[0 a 0];
      endo=1;
end
 
if frontsandelgs==1,
   % create wave fronts and resulting infinite medium potentials
   % at OBS and at Vera and VERb
   
   [VERAN,ITRIAN]=make_annulus(a,b,nr,nphi,0);
   VER=VERAN; ITRI=ITRIAN;
   nver=size(VER,1);
   figure(1)
   clf
   triplot
   
   for i=1:nver;
       DIST(i,1)=distance_annulus(pnt1,VER(i,:),a,b);
   end

   % ms if velocity in m/s
   TIMS=DIST*10/velocity;
   figure(1)
   clf
   cmap='tims.mcm';
   zebra=-0.25;
   grsw=0;
   lsw=0;
   iview=9;
   % triplot
   % hold on
   % identify isochrone: j
   nr=11;
   for i=1:nr,
      r(i)=b-(i-1)*(b-a)/(nr-1);
   end
   nphi=800;
   for i=1:nphi+1,,
      phi(i)=(i-1)*pi/nphi;
      x(i)=sin(phi(i));
      y(i)=cos(phi(i));
      for j=1:nr,
          D(i,j)=distance_annulus(pnt1,r(j)*[x(i)+eps y(i),0],a,b);
      end
   end
   
   %figure(3)
   %clf
   %plot(D)

   step=0.05;
   figure(1)
   d=step:step:max(max(D));
   nd=length(d);
   RES=zeros(nd,nr);
   X=RES;
   Y=RES;
   for id=1:nd,
       for j=1:nr,
           if d(id)>=D(1,j) & d(id) < D(nphi+1,j),
             is=find(D(:,j) <=d(id));
             if isempty(is)==0, 
                ir=reverse(is);
                RES(id,j)=phi(ir(1));
                X(id,j)=r(j)*cos(RES(id,j));
                Y(id,j)=r(j)*sin(RES(id,j));
             end
           end
       end
   end

   % treat incomplete number of nodes identified on an isochrone; 
   % force all to carry nr nodes

   for id=1:nd,
       if RES(id,1)~=0 | RES(id,nr)~=0,
          rho=d(id);
          %['refine' num2str(id)]
          if endo==1,
             % endocardial focus 
             if RES(id,1)==0,
                alpha=atan2(r(nr)*sin(RES(id,nr))-pnt1(1),r(nr)*cos(RES(id,nr))-pnt1(2));
                beta=(0:nr-2)/(nr-1)*alpha;
                XY=[ cos(beta)*rho; sin(beta)*rho;];
                X(id,1:nr-1)=XY(1,:)+ones(1,nr-1)*pnt1(2);
                Y(id,1:nr-1)=XY(2,:)+ones(1,nr-1)*pnt1(1);
             end
          else,
             % epicardial focus
             if RES(id,nr)==0 & id < nd/2,,
                alpha=atan2(r(1)*sin(RES(id,1))-pnt1(1),-r(1)*cos(RES(id,1))+pnt1(2));
                beta=(nr-1:-1:0)/(nr-1)*alpha;
                XY=[-cos(beta)*rho; sin(beta)*rho;];
                X(id,1:nr)=XY(1,:)+ones(1,nr)*pnt1(2);
                Y(id,1:nr)=XY(2,:)+ones(1,nr)*pnt1(1);
             end
          end
       
          if RES(id,1)~=0 | RES(id,nr)~=0,
             rho=d(id);
             if endo==1,
                % endocardial focus 
                if RES(id,1)==0,
                   alpha=atan2(r(nr)*sin(RES(id,nr))-pnt1(1),r(nr)*cos(RES(id,nr))-pnt1(2));
                   beta=(0:nr-2)/(nr-1)*alpha;
                   XY=[ cos(beta)*rho; sin(beta)*rho;];
                   X(id,1:nr-1)=XY(1,:)+ones(1,nr-1)*pnt1(2);
                   Y(id,1:nr-1)=XY(2,:)+ones(1,nr-1)*pnt1(1);
                end
             else,
                % epicardial focus
                if RES(id,nr)==0 & id < nd/2,,
                   alpha=atan2(r(1)*sin(RES(id,1))-pnt1(1),-r(1)*cos(RES(id,1))+pnt1(2));
                   beta=(nr-1:-1:0)/(nr-1)*alpha;
                   XY=[-cos(beta)*rho; sin(beta)*rho;];
                   X(id,1:nr)=XY(1,:)+ones(1,nr)*pnt1(2);
                   Y(id,1:nr)=XY(2,:)+ones(1,nr)*pnt1(1);
                end
             end
       
             if RES(id,nr)==0 & id > nd/2,
                list=find(RES(id,:)==0);
                nl=list(1)-1;
                %'refine arc'
                Y(id,nl+1)=0;
                if nl>=2,
                   aa=Y(id,nl-1)/Y(id,nl);
                   X(id,nl+1)=(X(id,nl-1)-aa*X(id,nl))/(1-aa);
                end
                if nl==1,
                   X(id,nl+1)=-b+2*(X(id,nl)+b);
                end
                x=X(id,1:nl+1);
                xx=x(1)+(0:nr-1)/(nr-1)*(x(nl+1)-x(1));
                y=spline(x,Y(id,1:nl+1),xx);
                X(id,:)=xx;
                Y(id,:)=y;             
             end
          end
       end
   end
 
   nfronts=size(X,1);
   for i=1:nfronts;
       plot3(-Y(i,:),X(i,:),zeros(1,nr),'W')
   end

   % create triangle specs wave front
   nphi=24;
   phi=(0:nphi-1)/nphi*2*pi;
   ntri=0;
   ITRIF=[];
   for j=1:nr-1;
       for i=1:nphi-1,
           ntri=ntri+1;
           ITRIF(ntri,1:3)=[ (j-1)*nphi+i  (j-1)*nphi+i+1 j*nphi+i];
           ntri=ntri+1;
           ITRIF(ntri,1:3)=[ (j-1)*nphi+i+1 j*nphi+i+1 j*nphi+i];
       end
       ntri=ntri+1;
       ITRIF(ntri,1:3)=[ j*nphi (j-1)*nphi+1  (j+1)*nphi];
       ntri=ntri+1;
       ITRIF(ntri,1:3)=[ (j-1)*nphi+1 j*nphi+1 (j+1)*nphi];
   end
   
   % create triangle specs cross-section myocard
   ntri=0;
   ITRIC=[];
   for j=1,
       for i=1:nphi-1,
           ntri=ntri+1;
           ITRIC(ntri,1:3)=[ (j-1)*nphi+i  (j-1)*nphi+i+1 j*nphi+i];
           ntri=ntri+1;
           ITRIC(ntri,1:3)=[ (j-1)*nphi+i+1 j*nphi+i+1 j*nphi+i];
       end
       ntri=ntri+1;
       ITRIC(ntri,1:3)=[ j*nphi (j-1)*nphi+1  (j+1)*nphi];
       ntri=ntri+1;
       ITRIC(ntri,1:3)=[ (j-1)*nphi+1 j*nphi+1 (j+1)*nphi];
   end
   
   figure(1)
   hold on
   iview=1;
   
   for t=1:nfronts, 
       clf
       nver=0;
       VERF=[];
       for j=nr:-1:1,
           rho=Y(t,j); z=X(t,j);
           for k=1:nphi,
               nver=nver+1;
               VERF(nver,1:3)=[-rho*cos(phi(k)) rho*sin(phi(k)) z];
           end
       end
    
       edgea=[1:nphi 1];
       edgeb=[nver-nphi+1:nver nver-nphi+1];
   
       % make boundary contour
       nphi=24;
       phi=(0:nphi-1)/nphi*2*pi;
       Bb=[zeros(nphi+1,1) rb*[cos(phi) cos(phi(1))]' rb*[sin(phi) sin(phi(1))]'];
       Ba=ra/rb*Bb;
       VERC=[Bb(1:nphi,:); Ba(1:nphi,:)];
       [VER,ITRI]=addtwotris(VERC,ITRIC,VERF,ITRIF);
       VALS=[];
       lsw=1;
       zebra=[];
       boxsize=4.5;
       triplot
       view(90,-60)
       hold on
       plot3(Bb(:,1),Bb(:,2),Bb(:,3),'b')
       plot3(Ba(:,1),Ba(:,2),Ba(:,3),'r')
       set(ht,'vis','off')
       
       next=1;
       if next==1,
           % used for switching off time markers
          uitim=uicontrol('style','text');
          uiboxtim=[.18 .92 .15 .05];
          set(uitim,'units','norm','position',uiboxtim',...
          'string',['t= ' num2str(t*10*step)],'fontsize',10)
       end
       
       % infinite medium potentials at OBS
       for kk=1:size(OBS,1),
           SA=dsa(VERF,ITRIF,OBS(kk,1:3),0.1);
           PHI(kk,t)=-10*sum(sum(SA))/pi;
       end
    
       % infinite medium potentials at vertices of endocard 
       for kk=1:size(VERa,1),
           SA=dsa(VERF,ITRIF,VERa(kk,1:3),0.1);
           PHIa(kk,t)=-10*sum(sum(SA))/pi;
       end
    
       % infinite medium potentials at vertices of epicard 
       for kk=1:size(VERa,1),
           SA=dsa(VERF,ITRIF,VERb(kk,1:3),0.1);
           PHIb(kk,t)=-10*sum(sum(SA))/pi;
       end
    
       plot3(OBS(1:nobs,1),OBS(1:nobs,2),OBS(1:nobs,3),'k+')
    
       next=0;
       if next==1,
          %used for illustrating field points
          delete(hcbar)
          plot3(OBS(1:nobs,1),OBS(1:nobs,3),OBS(1:nobs,2),'k+')
          text(0,3.5,2.2,'A3','fontsize',12)
          text(0,4,1,'A2','fontsize',12)
          text(0,1.7,0.65,'B2','fontsize',12)
          text(0,1.2,1.35,'B4','fontsize',12)
          text(0,3.1,-.25,'K','fontsize',12)
          text(0,0.2,3.2,'L','fontsize',12)
          text(0,-3.1,-0.25,'M','fontsize',12)
          text(0,0,-0.25,'C','fontsize',12)
          pause
       end
       pause(0.6)
   end
   
   nt=size(PHI,2);
   % 200=100 ms if v=1 m/s 
   nextra=210-nt;
   PHI=[zeros(nobs,10) PHI zeros(nobs,nextra)];
   PHIobs=PHI;
   nv=size(VERa,1);
   PHIa=[zeros(nv,10) PHIa zeros(nv,nextra)];
   nv=size(VERb,1);
   PHIb =[zeros(nv,10) PHIb zeros(nv,nextra)];
   if endo==1,
       savemat('phiobs_endostim.mat',PHIobs)
       savemat('phiendo_endostim.mat',PHIa)
       savemat('phiepi_endostim.mat',PHIb)
   else,
       savemat('phiobs_epistim.mat',PHIobs)
       savemat('phiendo_epistim.mat',PHIa)
       savemat('phiepi_epistim.mat',PHIb)
   end
   % end of preliminary computation of infinite medium potentials
end 
   
if endo==1,
   PHIobs_inf=loadmat('phiobs_endostim.mat');
   PHIa_inf  =loadmat('phiendo_endostim.mat');
  
else,
   PHIobs_inf=loadmat('phiobs_epistim.mat');
   PHIa_inf  =loadmat('phiendo_epistim.mat');
end

[VERa,ITRIa]=loadtri('endocard.tri');
nvera=length(VERa);
if ~exist('B'),
   for i=1:nvera,
       [B(i,1:nvera)]=rowforw(VERa,ITRIa,VERa(i,:));
   end
end

siga=3; % cavity conductivity (relative to myocard and exterior)
if siga~=sigp, inhom=1, end
if inhom==1,
    'inhom'
    % compute matrices for treating boundary effects 
    % below: Brody effect only;: boundary: VERa
  
    sigp=1;
    sigm=siga;
    sigs=1;
    
    kappa=(sigm-sigp)/(sigm+sigp);
    C=inv(eye(nvera)-kappa*B); % use to compute boundary potentials
    nobs=size(OBS,1);
    PHIa_inh=2*sigs/(sigp+sigm)*C*PHIa_inf;
    
    % compute the transfer from secondary sources at boundary to observation points
    locate=zeros(1,nobs);
    for i=1:nobs,
        [T(i,1:nvera),jsing(i)]=rowforw(VERa,ITRIa,OBS(i,:));
        locate(i)=sum(T(i,:),2); 
    end
   
% create the total solution at OBS
  nvera=length(VERa);
  
  for i=1:nobs,
      sigobs(i)=sigm; if locate(i) <.1, sigobs(i)=sigp; end
      PHIobs_inf_scaled(i,:)=sigs*PHIobs_inf(i,:)/sigobs(i);
      
      if locate(i)>1.5 | locate(i) < .5,
      % if jsing(i)==0,
         PHIobs_sec(i,:) = 0.5*(sigm-sigp)/sigobs(i)*T(i,:)*PHIa_inh;
         PHIobs_inh(i,:)=PHIobs_inf_scaled(i,:)+PHIobs_sec(i,:);  
       else,% obs on boundary; 
         PHIobs_inh(i,:)=PHIa_inh(jsing(i),:);
         PHIobs_sec(i,:)=PHIobs_inh(i,:)-PHIobs_inf_scaled(i,:);
      end
  end   
end

%[(1:nobs)' OBS(:,3) jsing' locate' sigobs']

t=(1:size(PHIa_inf,2))/2;

figure(1)
clf

iver=1;
plot(t,PHIa_inf(iver,:))
xlabel('time; (ms)')

hold on
plot([t(1) t(length(t))], [0 0], 'k:')
title([ 'potential at boundary vertex: ' num2str(iver)])
if exist('PHIa_inh'),plot(t,PHIa_inh(iver,:),'r'); end

imark=find(jsing~=0);
jmark=jsing(imark);

for fig=2:3,
    figure(fig)
    
    clf
    if fig==2, 
       showtim=8; % time from stimulus
    else,
       showtim=24;
    end
    tt=2*showtim+10;
    obs=1:length(OBS);
    plot(OBS(obs,3),PHIobs_inf(obs,tt),'-')
    hold on
    plot(OBS(imark,3),PHIobs_inf(obs(imark),tt),'k*')
    xlabel('distance along z-axis; (cm)')
    if endo ==1, 
       title(['time from stimulus = ' num2str((tt-10)/2) ' ms;  endo stim']),
    else,
       title(['time from stimulus = ' num2str((tt-10)/2) ' ms;  epi stim']),
    end
    plot([OBS(1,3) OBS(length(OBS),3)],[0 0],'k:')

    if exist('PHIobs_inh'),
       plot(OBS(obs,3),   PHIobs_inh(obs,tt),'r');
       plot(OBS(imark,3), PHIobs_inh(obs(imark),tt),'k*')
    
    plot(OBS(obs,3),     PHIobs_sec(obs,tt),'g'); 
    plot(OBS(imark,3),   PHIobs_sec(obs(imark),tt),'k*')
    
    plot(OBS(obs,3),   PHIobs_inf_scaled(obs,tt),'m'); 
    plot(OBS(imark,3), PHIobs_inf_scaled(obs(imark),tt),'k*')
    plot(OBS(imark,3), PHIa_inh(jmark,tt),'r*')
    
    end
end

figure(4)
clf
ipl=0;
for showtime=6:4:24,
    ipl=ipl+40;
    tt=2*showtime+1;
    obs=1:length(OBS);
    plot(OBS(obs,3),ipl+PHIobs_inf(obs,tt),'b:')
    hold on
    %plot(OBS(imark,3),ipl+PHIobs_inf(obs(imark),tt),'k*')
    set(gcf,'color',[1 1 1])
    %axis off
    axis([0 4.5 0 220])
    box off
    xlabel('distance along line of IME; (cm)')
    ylabel('potential; profiles spaced at 40 mV')
    if endo ==1, 
       title('potential profile; endo stimulus'),
    else,
        title('potential profile; epi stimulus'),
    end
    plot([OBS(1,3) OBS(length(OBS),3)],ipl+[0 0],'k:')
    plot([2.4 4],ipl+[0 0],'k', 'linewidth',1.5)
    text(4.5,ipl-6,num2str((tt-9)/2))
    if exist('PHIobs_inh'),
       plot(OBS(obs,3),   ipl+PHIobs_inh(obs,tt),'r');
       plot(OBS(imark,3), ipl+PHIobs_inh(obs(imark),tt),'k+')
    end
end
     plot([2.4 4],[0 0],'k', 'linewidth',1.5)
     axis off
     
     plot([0 4.5],[0 0],'k')
     % create and ticks
     for i=1:5,
         plot([i-1 i-1],[0 3],'k')
         text(i-1.025,-8,num2str(i-1))
     end
     

