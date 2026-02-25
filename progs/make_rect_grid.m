% make_rect_grid.m
% function [VER, ITRI]=make_rect grid(nx,ny);
% generates mesh for unit squares with edges from origin to 
% [1 0 0] [0 1 0] with nodes nx and ny, respectively
% hence: its normal is oriented along z-axis

% 2017_09_27; A. van Oosterom


n_y=11;  
n_x=8;

VER=[];
for j=1:n_x+1,
    for i=1:n_y+1,
        VER=[VER;[j i]];
    end
end

SQUARE=[1 n_y+2 n_y+3 2];
for i=2:n_y,
    SQUARE=[SQUARE; SQUARE(end,:) + ones(1,4)];
end

for i=2:n_x,
    SQUARE=[SQUARE;[SQUARE(end-n_y+1:end,:)+ (n_y+1)*ones(n_y,4)]];
end

n_squares=n_y*n_x;

figure(1)
clf
hold on  

for i=1:n_y+1;
    plot([1 n_x+1], [i i]);
    hold on   
end

hold on  
for i=1:n_x+1;
    plot([i i],[1 n_y+1]);
    hold on   
end

axis off
axis square
axis equal


pause

for i=2:n_y-1;
    cent=mean(VER(SQUARE(i,:),:));
    plot(cent(1),cent(2),'*r')
    hold on   
end

for i=2:n_y-1;
    cent=mean(VER(SQUARE(i+n_y,:),:));
    plot(cent(1),cent(2),'*r')
    hold on   
end








