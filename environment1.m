%% environment1.m  --  Simple Rectangular Obstacle Map
% -------------------------------------------------------------------------
% 1) Domain and grid definition
% -------------------------------------------------------------------------

nx = 300;                         % grid points (x)
ny = 300;                         % grid points (y)
xrange = [0, 10];                 % x‑axis limits
yrange = [0, 10];                 % y‑axis limits

xvec = linspace(xrange(1), xrange(2), nx);
yvec = linspace(yrange(1), yrange(2), ny);
[XX, YY] = meshgrid(xvec, yvec);  % XX(i,j)=x , YY(i,j)=y

obsMap = zeros(ny, nx);           % 0 = free, 1 = obstacle

% -------------------------------------------------------------------------
% 2) Define rectangular obstacle  [xmin, ymin, width, height]
% -------------------------------------------------------------------------
rect = [3.0, 3.0, 7.0, 4.0];      % large horizontal block

mask = (XX >= rect(1) & XX <= rect(1)+rect(3)) & ...
       (YY >= rect(2) & YY <= rect(2)+rect(4));
obsMap(mask) = 1;

% -------------------------------------------------------------------------
% 3) Start / Goal positions
% -------------------------------------------------------------------------
startPos = [0, 0];                % bottom‑left
goalPos  = [10, 10];              % top‑right

% -------------------------------------------------------------------------
% 4) Visualization
% -------------------------------------------------------------------------
figure('Color','w');
imagesc(xvec, yvec, obsMap);
colormap([1 1 1; 0 0 0]); caxis([0 1]);
set(gca,'YDir','normal'); axis equal tight; hold on;
title('Environment  –  Rectangular Obstacle');

plot(startPos(1), startPos(2), 'bs', 'MarkerFaceColor','b', 'MarkerSize',8);
plot(goalPos(1),  goalPos(2),  'gs', 'MarkerFaceColor','r', 'MarkerSize',8);

text(startPos(1)+0.3, startPos(2)+0.3, 'Start', 'Color','b', 'FontWeight','bold');
text(goalPos(1)-1.2, goalPos(2)-0.4, 'Goal',  'Color','r', 'FontWeight','bold');

xlim(xrange); ylim(yrange);
