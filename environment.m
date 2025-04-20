%% 1) Domain and grid definition
clear; clc; close all;

nx = 300;   % number of grid points along x
ny = 300;   % number of grid points along y
xrange = [0, 10];  % x‑axis limits
yrange = [0, 10];  % y‑axis limits

xvec = linspace(xrange(1), xrange(2), nx);
yvec = linspace(yrange(1), yrange(2), ny);
[XX, YY] = meshgrid(xvec, yvec);   % XX(i,j)=x,  YY(i,j)=y

%% 2) Define circular obstacles:  [center‑x, center‑y, radius]
obstacles = [
    2, 3.5, 2.8;
    7, 1, 1.7;
    8, 6.5, 3.0;
    2, 9, 2.0
];

% Binary obstacle map: 1 = obstacle, 0 = free
obsMap = zeros(ny, nx);

for k = 1:size(obstacles,1)
    cx = obstacles(k,1);
    cy = obstacles(k,2);
    r  = obstacles(k,3);

    distSq = (XX - cx).^2 + (YY - cy).^2;
    inside = (distSq <= r^2);
    obsMap(inside) = 1;
end

%% 3) Start / goal positions
startPos = [0, 0];
goalPos  = [10, 10];

%% 4) Visualization
figure('Color','w');

% Custom colormap:
%   0 → background (white)
%   1 → obstacle   (black)
myColorMap = [1 1 1;   % background
              0 0 0];  % obstacle
colormap(myColorMap);

% Keep colorbar range in [0,1] so that the colormap indices line up
imagesc(xvec, yvec, obsMap);
caxis([0 1]);
set(gca,'YDir','normal');
axis equal tight; hold on;
title('Obstacle Map (black = obstacle, white = free space)');

% Plot start (blue) and goal (green with red fill) markers
plot(startPos(1), startPos(2), 'bo', 'MarkerFaceColor','b', 'MarkerSize',8);
plot(goalPos(1),  goalPos(2),  'go', 'MarkerFaceColor','r', 'MarkerSize',8);

% Annotate the start and goal without overlap
textOffset = 0.3;
text(startPos(1)+textOffset, startPos(2)+textOffset, 'Start', ...
     'Color','blue','FontWeight','bold', ...
     'BackgroundColor','w','Margin',2);
text(goalPos(1)-1.5, goalPos(2)-textOffset, 'Goal', ...
     'Color','red','FontWeight','bold', ...
     'BackgroundColor','w','Margin',2);

xlim(xrange); ylim(yrange);
