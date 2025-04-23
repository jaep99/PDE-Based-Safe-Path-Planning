function QualitativeAnalysis()
% QualitativeAnalysis: Shows paths for different (alpha, beta) values 
% and environments in a 3x2 tiled layout.
%
% - Closes existing figures and sets a default window size.
% - Iterates over 3 environments and 2 parameter pairs.
% - For each combination, computes a path using runPDE_OneShot 
%   and then plots the obstacle map, start/goal, and path.
% - Uses a tiled layout with compact spacing.

    % Close figures, set default size/position
    close all;
    defaultPos = [100, 100, 1200, 1000];
    set(0, 'DefaultFigurePosition', defaultPos, 'DefaultFigureUnits', 'pixels');

    % Create main figure and center it
    fig = figure('Name', 'Qualitative combos', 'NumberTitle', 'off', 'Color', 'w');
    movegui(fig, 'center');

    % Tiled layout for 3 environments x 2 parameter sets
    layout = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Define environment scripts and (alpha,beta) pairs
    envList = {'environment1', 'environment2', 'environment3'};
    alphaBetaPairs = [
        5, 0.3;  % alpha=5, beta=0.3
        2, 1.0;  % alpha=2, beta=1.0
    ];

    % Loop through environments
    for eIdx = 1:numel(envList)
        % Load environment variables (obsMap, xvec, yvec, startPos, goalPos)
        run(envList{eIdx});

        % Close any extra figure that might have been opened
        figs = findobj('Type', 'figure');
        for k = 1:numel(figs)
            if figs(k) ~= fig
                close(figs(k));
            end
        end

        % Loop through each (alpha, beta) combo
        for cIdx = 1:size(alphaBetaPairs, 1)
            alpha = alphaBetaPairs(cIdx, 1);
            beta  = alphaBetaPairs(cIdx, 2);

            % Compute path and metrics
            [pathXY, finalCost, pathLen] = runPDE_OneShot(obsMap, startPos, ...
                goalPos, xvec, yvec, alpha, beta);

            % Move to next tile
            ax = nexttile(layout);

            % Display obstacle map
            imagesc(xvec, yvec, obsMap, 'Parent', ax);
            colormap(ax, [1 1 1; 0 0 0]);
            set(ax, 'YDir', 'normal', 'DataAspectRatio', [1 1 1]);
            xlim(ax, [min(xvec), max(xvec)]);
            ylim(ax, [min(yvec), max(yvec)]);
            hold(ax, 'on');

            % Plot start, goal, and path
            plot(ax, startPos(1), startPos(2), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
            plot(ax, goalPos(1),  goalPos(2),  'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
            plot(ax, pathXY(:,1), pathXY(:,2), 'r-', 'LineWidth', 2);

            % Title
            greekAlpha = char(945);  % α
            greekBeta  = char(946);  % β
            title(ax, sprintf('%s=%.1f, %s=%.1f | Cost=%.1f, Len=%.1f', ...
                greekAlpha, alpha, greekBeta, beta, finalCost, pathLen));
            hold(ax, 'off');
        end
    end
end


% -- Helper function --
function [pathXY, finalCost, pathLen] = runPDE_OneShot(obsMap, startPos, goalPos, xvec, yvec, alpha, beta)
% runPDE_OneShot: Single-run Eikonal PDE solver for given alpha, beta.
% Returns path (pathXY), cost (finalCost), and path length (pathLen).

    distToObs = bwdist(logical(obsMap));
    cMap = 1 + alpha * exp(-beta * distToObs);
    cMap(obsMap == 1) = 1e5;

    [ny, nx] = size(obsMap);
    T = inf(ny, nx);

    [sI, sJ] = xy2ij(startPos, xvec, yvec);
    [gI, gJ] = xy2ij(goalPos,  xvec, yvec);
    T(sI, sJ) = 0;

    state = zeros(ny, nx, 'uint8');
    state(sI, sJ) = 1;

    pq.data = zeros(0, 2);
    pq.keys = [];
    pq = pushOrUpdatePQ(pq, [sI, sJ], 0);

    neighbors = [
       -1  0;  1  0;  0 -1;  0  1;
       -1 -1; -1  1;  1 -1;  1  1
    ];

    % Fast-marching
    while ~isempty(pq.data)
        [curr, val, pq] = popPQ(pq);
        ci = curr(1); cj = curr(2);
        if state(ci,cj) == 2, continue; end
        state(ci,cj) = 2;
        if ci == gI && cj == gJ, break; end
        for nn = 1:size(neighbors,1)
            ni = ci + neighbors(nn,1);
            nj = cj + neighbors(nn,2);
            if ni<1 || ni>ny || nj<1 || nj>nx, continue; end
            if state(ni,nj) == 2, continue; end
            cand = val + norm(neighbors(nn,:)) * cMap(ni,nj);
            if cand < T(ni,nj)
                T(ni,nj) = cand;
                pq = pushOrUpdatePQ(pq, [ni, nj], cand);
                state(ni,nj) = 1;
            end
        end
    end

    % Backtrace
    pathIJ = [gI, gJ];
    cur = [gI, gJ];
    while any(cur ~= [sI, sJ])
        ci = cur(1); cj = cur(2);
        bestPt = cur; bestVal = T(ci, cj);
        for nn = 1:size(neighbors,1)
            ni = ci + neighbors(nn,1);
            nj = cj + neighbors(nn,2);
            if ni<1 || ni>ny || nj<1 || nj>nx, continue; end
            if T(ni,nj) < bestVal
                bestVal = T(ni,nj);
                bestPt = [ni, nj];
            end
        end
        if all(bestPt == cur), break; end
        pathIJ = [pathIJ; bestPt];
        cur = bestPt;
    end

    pathXY = ij2xy(pathIJ, xvec, yvec);
    finalCost = T(gI, gJ);

    d = diff(pathXY, 1, 1);
    pathLen = sum(sqrt(d(:,1).^2 + d(:,2).^2));
end

function pq = pushOrUpdatePQ(pq, item, key)
% pushOrUpdatePQ: Inserts or updates an item in the priority queue.
    idx = find(ismember(pq.data, item, 'rows'));
    if ~isempty(idx)
        if key < pq.keys(idx), pq.keys(idx) = key; end
    else
        pq.data = [pq.data; item]; 
        pq.keys = [pq.keys; key];
    end
end

function [item, val, pq] = popPQ(pq)
% popPQ: Removes and returns the item with the smallest key.
    [~, k] = min(pq.keys);
    item = pq.data(k,:);
    val  = pq.keys(k);
    pq.data(k,:) = [];
    pq.keys(k)   = [];
end

function [i, j] = xy2ij(xy, xvec, yvec)
% xy2ij: Converts (x,y) to nearest grid indices.
    [~, j] = min(abs(xvec - xy(1)));
    [~, i] = min(abs(yvec - xy(2)));
end

function xy = ij2xy(ijlist, xv, yv)
% ij2xy: Converts (i,j) to (x,y) coordinates.
    xy = nan(size(ijlist,1), 2);
    for n = 1:size(ijlist,1)
        xy(n,:) = [xv(ijlist(n,2)), yv(ijlist(n,1))];
    end
end
