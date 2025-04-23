function [finalCost, pathLen, runtime] = EikonalPDE_result()
% Solves the Eikonal PDE ||grad T|| = c(x) with a soft cost near obstacles.
% Then it backtraces from the goal to the start to obtain a global path.
%
% Steps:
%   1) Check environment variables (obsMap, startPos, goalPos, xvec, yvec, nx, ny).
%   2) Define cMap = 1 + alpha * exp(-beta * distToObs), obstacles = 1e5.
%   3) Solve T with a simple fast-marching approach (starting from T=0 at start).
%   4) Backtrace from the goal to the start to form a path.
%   5) Plot the final path.
%
% Outputs:
%   finalCost : T(gI, gJ), minimal cost from start to goal
%   pathLen   : Euclidean length of the path
%   runtime   : total runtime in seconds for PDE + backtrace

tic; % start timer

%% 1) Check environment
whosEnv = evalin('base','whos');
envVars = {whosEnv.name};
reqVars = {'obsMap','startPos','goalPos','xvec','yvec','nx','ny'};
for rv = reqVars
    if ~ismember(rv{1}, envVars)
        error(['Variable ', rv{1}, ' not found in base workspace. Please run environment script first.']);
    end
end

obsMap   = evalin('base','obsMap');
startPos = evalin('base','startPos');
goalPos  = evalin('base','goalPos');
xvec     = evalin('base','xvec');
yvec     = evalin('base','yvec');
nx       = evalin('base','nx');
ny       = evalin('base','ny');

%% 2) Build soft cost map
distToObs = bwdist(logical(obsMap));
alpha = 5;
beta  = 0.3;
cMap = 1 + alpha * exp(-beta * distToObs);
cMap(obsMap == 1) = 1e5; % large cost for obstacles

%% 3) Fast-marching solution for Eikonal
T = inf(ny, nx);
[sI, sJ] = xy2ij(startPos, xvec, yvec);
[gI, gJ] = xy2ij(goalPos,  xvec, yvec);
T(sI, sJ) = 0;

state = zeros(ny, nx, 'uint8'); % 0=far, 1=trial, 2=accepted
state(sI, sJ) = 1;

% Priority queue
pq.data = zeros(0, 2);
pq.keys = [];

% Push start node
pq = pushOrUpdatePQ(pq, [sI, sJ], 0);

neighbors = [
   -1  0;
    1  0;
    0 -1;
    0  1;
   -1 -1;
   -1  1;
    1 -1;
    1  1
];

while ~isempty(pq.data)
    [curr, val, pq] = popPQ(pq);
    ci = curr(1); cj = curr(2);
    if state(ci,cj) == 2, continue; end
    state(ci,cj) = 2; % accepted

    if ci == gI && cj == gJ
        break; % goal found
    end

    for nn = 1:size(neighbors,1)
        ni = ci + neighbors(nn,1);
        nj = cj + neighbors(nn,2);
        if ni<1 || ni>ny || nj<1 || nj>nx, continue; end
        if state(ni,nj) == 2, continue; end

        cost = cMap(ni,nj);
        stepd = norm(neighbors(nn,:)); % 1 or sqrt(2)
        cand = val + stepd * cost;
        if cand < T(ni,nj)
            T(ni,nj) = cand;
            pq = pushOrUpdatePQ(pq, [ni, nj], cand);
            state(ni,nj) = 1;
        end
    end
end

%% 4) Backtrace path
pathIJ = [gI, gJ];
cur = [gI, gJ];
while any(cur ~= [sI, sJ])
    ci = cur(1); cj = cur(2);
    bestPt = cur; bestVal = T(ci, cj);
    for nn = 1:size(neighbors,1)
        ni = ci + neighbors(nn,1);
        nj = cj + neighbors(nn,2);
        if ni<1 || ni>ny || nj<1 || nj>nx, continue; end
        tv = T(ni,nj);
        if tv < bestVal
            bestVal = tv;
            bestPt = [ni, nj];
        end
    end
    if all(bestPt == cur)
        break; % stuck
    end
    pathIJ = [pathIJ; bestPt];
    cur = bestPt;
end

pathXY = ij2xy(pathIJ, xvec, yvec);

%% Compute final cost and path length
finalCost = T(gI, gJ);
diffs = diff(pathXY, 1, 1);
segLens = sqrt(diffs(:,1).^2 + diffs(:,2).^2);
pathLen = sum(segLens);

runtime = toc; % end timer

%% Print results
fprintf('Results:\n');
fprintf('  alpha=%.2f, beta=%.2f\n', alpha, beta);
fprintf('  Final cost:    %.4f\n', finalCost);
fprintf('  Path length:   %.4f\n', pathLen);
fprintf('  Runtime (sec): %.4f\n\n', runtime);

%% 5) Plot the path
figure('Name','Eikonal Path','Color','w');
imagesc(xvec, yvec, obsMap);
colormap([1 1 1; 0 0 0]);
set(gca, 'YDir', 'normal');
axis equal tight; hold on;
plot(startPos(1), startPos(2), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot(goalPos(1), goalPos(2), 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(pathXY(:,1), pathXY(:,2), 'r-', 'LineWidth', 2);
title('Eikonal-based Path');
hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%
function pq = pushOrUpdatePQ(pq, item, key)
% pushOrUpdatePQ: Inserts or updates an item in the priority queue.
idx = find(ismember(pq.data, item, 'rows'));
if ~isempty(idx)
    if key < pq.keys(idx)
        pq.keys(idx) = key;
    end
else
    pq.data = [pq.data; item];
    pq.keys = [pq.keys; key];
end
end

function [item, val, pq] = popPQ(pq)
% popPQ: Removes and returns the item with the smallest key.
if isempty(pq.data)
    item = [];
    val  = inf;
    return;
end
[~, k] = min(pq.keys);
item = pq.data(k,:);
val  = pq.keys(k);
pq.data(k,:) = [];
pq.keys(k)   = [];
end

function [i, j] = xy2ij(xy, xvec, yvec)
% xy2ij: Converts (x,y) to (i,j) based on the nearest grid indices.
x = xy(1); y = xy(2);
[~, j] = min(abs(xvec - x));
[~, i] = min(abs(yvec - y));
end

function xy = ij2xy(ijlist, xv, yv)
% ij2xy: Converts a list of (i,j) indices to (x,y) positions.
xy = zeros(size(ijlist,1), 2);
for n = 1:size(ijlist,1)
    i = ijlist(n,1);
    j = ijlist(n,2);
    if i<1 || j<1 || i>numel(yv) || j>numel(xv)
        xy(n,:) = [NaN, NaN];
    else
        xy(n,:) = [xv(j), yv(i)];
    end
end
end
