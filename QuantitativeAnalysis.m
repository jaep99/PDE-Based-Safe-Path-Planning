function QuantitativeAnalysis()
% QuantitativeAnalysis: Performs a parameter sweep of (alpha, beta)
% over multiple environments, then prints [Env, Alpha, Beta, Cost, Length, Runtime].
%
% - Defines environment list (environment1, environment2, environment3).
% - Defines alphaList, betaList.
% - For each combination, runs the Eikonal PDE solver without plotting.
% - Prints the results in a table-like format.

    envList   = {'environment1', 'environment2', 'environment3'};
    alphaList = [2, 5, 7, 9];
    betaList  = [0.3, 0.5, 0.7, 0.9];

    nEnv   = numel(envList);
    nAlpha = numel(alphaList);
    nBeta  = numel(betaList);
    total  = nEnv * nAlpha * nBeta;

    % Prepare storage for results (optional usage)
    results = zeros(total, 6);
    idx = 0;

    % Print header
    fprintf('Env\tAlpha\tBeta\tCost\tLength\tRuntime\n');

    % Sweep over environments and parameters
    for e = 1:nEnv
        % Load environment variables: obsMap, startPos, goalPos, xvec, yvec
        run(envList{e});

        for ai = 1:nAlpha
            for bi = 1:nBeta
                idx = idx + 1;
                a = alphaList(ai);
                b = betaList(bi);

                % Solve PDE
                [cost, len, t] = runEikonalNoPlot(obsMap, startPos, ...
                    goalPos, xvec, yvec, a, b);

                % Store and print
                results(idx, :) = [e, a, b, cost, len, t];
                fprintf('%d\t%.1f\t%.1f\t%.2f\t%.2f\t%.4f\n', e, a, b, cost, len, t);
            end
        end
    end
end


function [finalCost, pathLen, runtime] = runEikonalNoPlot(obsMap, startPos, goalPos, xvec, yvec, alpha, beta)
% runEikonalNoPlot: Solves Eikonal PDE without plotting; returns cost, path length, runtime.

    tic; % start timing

    % Build cost map
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

    pq.data = zeros(0,2);
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
        if state(ci, cj) == 2, continue; end
        state(ci, cj) = 2;
        if ci == gI && cj == gJ, break; end
        for k = 1:size(neighbors,1)
            ni = ci + neighbors(k,1);
            nj = cj + neighbors(k,2);
            if ni<1 || ni>ny || nj<1 || nj>nx, continue; end
            if state(ni,nj) == 2, continue; end
            cand = val + norm(neighbors(k,:)) * cMap(ni,nj);
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
        for k = 1:size(neighbors,1)
            ni = ci + neighbors(k,1);
            nj = cj + neighbors(k,2);
            if ni<1 || ni>ny || nj<1 || nj>nx, continue; end
            if T(ni, nj) < bestVal
                bestVal = T(ni, nj);
                bestPt = [ni, nj];
            end
        end
        if all(bestPt == cur), break; end
        pathIJ = [pathIJ; bestPt];
        cur = bestPt;
    end

    % Convert indices to XY and compute metrics
    pathXY = ij2xy(pathIJ, xvec, yvec);
    finalCost = T(gI, gJ);
    d = diff(pathXY, 1, 1);
    pathLen = sum(sqrt(d(:,1).^2 + d(:,2).^2));
    runtime = toc; % end timing
end

% Helper functions
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
    [~, k] = min(pq.keys);
    item = pq.data(k,:);
    val  = pq.keys(k);
    pq.data(k,:) = [];
    pq.keys(k)   = [];
end

function [i, j] = xy2ij(xy, xvec, yvec)
% xy2ij: Converts (x,y) to the closest (i,j) grid indices.
    [~, j] = min(abs(xvec - xy(1)));
    [~, i] = min(abs(yvec - xy(2)));
end

function xy = ij2xy(ijlist, xv, yv)
% ij2xy: Converts (i,j) grid indices to (x,y) coordinates.
    xy = nan(size(ijlist,1),2);
    for n = 1:size(ijlist,1)
        xy(n,:) = [xv(ijlist(n,2)), yv(ijlist(n,1))];
    end
end
