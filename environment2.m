%% environment2.m  --  Multiple Circular Obstacles Map
% -------------------------------------------------------------------------
% 1) Domain and grid definition
% -------------------------------------------------------------------------

nx = 300; ny = 300;
xrange = [0, 10]; yrange = [0, 10];

xvec = linspace(xrange(1), xrange(2), nx);
yvec = linspace(yrange(1), yrange(2), ny);
[XX, YY] = meshgrid(xvec, yvec);

obsMap = zeros(ny, nx);           % 0 = free, 1 = obstacle

% -------------------------------------------------------------------------
% 2) Define circular obstacles  [cx, cy, r]
% -------------------------------------------------------------------------
obstacles = [
    2, 3.5, 2.8;
    7, 1.0, 1.7;
    8, 6.5, 3.0;
    2, 9.0, 2.0
];

for k = 1:size(obstacles,1)
    cx = obstacles(k,1); cy = obstacles(k,2); r = obstacles(k,3);
    mask = (XX-cx).^2 + (YY-cy).^2 <= r^2;
    obsMap(mask) = 1;
end

% -------------------------------------------------------------------------
% 3) Start / Goal positions
% -------------------------------------------------------------------------
startPos = [0, 0];
goalPos  = [10, 10];

% -------------------------------------------------------------------------
% 4) Visualization
% -------------------------------------------------------------------------
figure('Color','w');
imagesc(xvec, yvec, obsMap);
colormap([1 1 1; 0 0 0]); caxis([0 1]);
set(gca,'YDir','normal'); axis equal tight; hold on;
title('Environment  –  Circular Obstacles');

plot(startPos(1), startPos(2), 'bs', 'MarkerFaceColor','b', 'MarkerSize',8);
plot(goalPos(1),  goalPos(2),  'gs', 'MarkerFaceColor','r', 'MarkerSize',8);

text(startPos(1)+0.3, startPos(2)+0.3, 'Start', 'Color','b', 'FontWeight','bold');
text(goalPos(1)-1.2, goalPos(2)-0.4, 'Goal',  'Color','r', 'FontWeight','bold');

xlim(xrange); ylim(yrange);
