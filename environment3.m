%% environment3.m  --  Maze Map
% -------------------------------------------------------------------------
% 1) Domain and grid definition
% -------------------------------------------------------------------------
rng(42);

nx = 300;                      % grid points (x)
ny = 300;                      % grid points (y)
xrange = [0, 10];              % x-axis limits
yrange = [0, 10];              % y-axis limits

xvec = linspace(xrange(1), xrange(2), nx);
yvec = linspace(yrange(1), yrange(2), ny);
[XX, YY] = meshgrid(xvec, yvec);

obsMap = zeros(ny, nx);        % 0 = free, 1 = obstacle
dx = xvec(2) - xvec(1);        % grid spacing

% -------------------------------------------------------------------------
% 2) Maze generation
% -------------------------------------------------------------------------
% Maze parameters
nCell = 17;                    % odd number of cells
cellM = 0.6;                   % cell size [m]
w     = 0.15;                  % wall thickness [m]
pW    = round(w / (2*dx));     % half wall (pixels)
cellP = round(cellM / dx);     % cell pitch (pixels)

% Helper function: clip indices to 1..nx or 1..ny
clip = @(idx, maxv) max(1, min(idx, maxv));

% Start with a full wall grid
obsMap(:) = 1;

% Compute the maze starting point (centered)
sx = clip(round((nx - (nCell-1)*cellP) / 2), nx);
sy = clip(round((ny - (nCell-1)*cellP) / 2), ny);

% Carve out cell centers
for i = 0:nCell-1
    for j = 0:nCell-1
        cx = clip(sx + i*cellP, nx);
        cy = clip(sy + j*cellP, ny);
        obsMap(clip(cy-pW, ny):clip(cy+pW, ny), ...
               clip(cx-pW, nx):clip(cx+pW, nx)) = 0;
    end
end

% Depth-first search (DFS) to carve out corridors
dxy = [1 0; -1 0; 0 1; 0 -1];  % E/W/N/S
visited = false(nCell, nCell);
stack = [1, 1];

while ~isempty(stack)
    curr = stack(end, :);
    visited(curr(2), curr(1)) = true;
    nbr = [];

    % Find unvisited neighbors
    for k = 1:4
        nxt = curr + dxy(k, :);
        if all(nxt >= 1 & nxt <= nCell) && ~visited(nxt(2), nxt(1))
            nbr = [nbr; k, nxt]; %#ok<AGROW>
        end
    end

    % Backtrack if no unvisited neighbors
    if isempty(nbr)
        stack(end, :) = [];
        continue;
    end

    % Pick a random neighbor
    pick = nbr(randi(size(nbr, 1)), :);
    dirk = pick(1); 
    nxt  = pick(2:3);

    % Convert cell coords to pixel coords
    cx = clip(sx + (curr(1)-1)*cellP, nx);
    cy = clip(sy + (curr(2)-1)*cellP, ny);
    nx_cx = clip(sx + (nxt(1)-1)*cellP, nx);
    nx_cy = clip(sy + (nxt(2)-1)*cellP, ny);

    % Carve out the passage
    if dirk <= 2  % E/W
        xRange = sort([cx, nx_cx]);
        obsMap(clip(cy-pW, ny):clip(cy+pW, ny), ...
               clip(xRange(1), nx):clip(xRange(2), nx)) = 0;
    else          % N/S
        yRange = sort([cy, nx_cy]);
        obsMap(clip(yRange(1), ny):clip(yRange(2), ny), ...
               clip(cx-pW, nx):clip(cx+pW, nx)) = 0;
    end

    stack(end+1, :) = nxt;
end

% -------------------------------------------------------------------------
% 3) Start / Goal positions
% -------------------------------------------------------------------------
startPos = [xvec(clip(sx, nx)),                    ...
            yvec(clip(sy, ny))];
goalPos  = [xvec(clip(sx + (nCell-1)*cellP, nx)),   ...
            yvec(clip(sy + (nCell-1)*cellP, ny))];

% -------------------------------------------------------------------------
% 4) Visualization
% -------------------------------------------------------------------------
figure('Color','w');
imagesc(xvec, yvec, obsMap);
colormap([1 1 1; 0 0 0]); 
caxis([0 1]);
set(gca, 'YDir', 'normal'); 
axis equal tight; hold on;
title(sprintf('Environment  –  DFS Maze', nCell, nCell));

plot(startPos(1), startPos(2), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot(goalPos(1),  goalPos(2),  'gs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

text(startPos(1)+0.3, startPos(2)+0.3, 'Start', 'Color','b', 'FontWeight','bold');
text(goalPos(1)-1.2, goalPos(2)-0.4, 'Goal',  'Color','r', 'FontWeight','bold');

xlim(xrange); ylim(yrange);
