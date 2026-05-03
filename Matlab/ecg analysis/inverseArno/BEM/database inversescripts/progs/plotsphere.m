function h = plotsphere(VER,ITRI,data)
   
   nver=size(VER,1);
   hs=patch('Vertices',VER,'Faces',ITRI,...
   'FaceVertexCData',ones(nver,3)*.7,'FaceColor','inter');
    set(hs,'EdgeColor',ones(1,3)*0.7,'Marker','.','MarkerEdgeColor','none');
    set(hs,'FaceAlpha',.7,'EdgeAlpha',.7);

axis vis3d

view(135,30)
hold on
fpi=2*pi/12;

% draw meridians  
for i=[1 2 4 5 7 8 10 11],
    p=[cos(i*fpi) sin(i*fpi) 0];
    plot_circle(p,50,0.5,'-','k');
    p=[0 0 cos(i*fpi-pi)];
    plot_circle(p,50,0.5,'-','k');
end

    equator=plot_circle([0 0 1.],50,1.5,'-','k');
    gmt=    plot_circle([0 1. 0],50,1.5,'-','k');
    left=   plot_circle([1. 0 0],50,1.5,'-','k');
    
    plot3([0;  1.], [0; 0], [0; 0],'k:','linewidth',1.5)
    plot3([1; 1.3], [0; 0], [0; 0],'k','linewidth',2)
    text(1.35,0,0,'front')
    
    plot3([0;  -1.], [0; 0], [0; 0],'k:','linewidth',1.5)
    plot3([-1; -1.3], [0; 0], [0; 0],'k','linewidth',2)
    text(-1.35,0,0,'back')
    
    plot3([0; 0], [0;   1], [0; 0],'k:','linewidth',1.5)
    plot3([0; 0], [1; 1.3], [0; 0],'k','linewidth',2)
    text(0,1.35,0,'left')
    plot3([0; 0], [0;   -1], [0; 0],'k:','linewidth',1.5)
    plot3([0; 0], [-1; -1.3], [0; 0],'k','linewidth',2)
    text(0,-1.35,0,'right')
    
    plot3([0; 0], [0; 0], [0;   1],'k:','linewidth',1.5)
    plot3([0; 0], [0; 0], [1; 1.3],'k','linewidth',2)
    text(0,0,1.35,'head')
    
    plot3([0; 0], [0; 0], [0;   -1],'k:','linewidth',1.5)
    plot3([0; 0], [0; 0], [1; -1.3],'k','linewidth',2)
    text(0,0,-1.35,'feet')
    
set(gca,'Visible','off')

%hold on

h = plot3(data(:,1),data(:,2),data(:,3),'r.');
set(h,'MarkerSize',14)