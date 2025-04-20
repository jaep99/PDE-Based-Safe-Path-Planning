function SnakeActiveContourImplementation()
% SnakeActiveContourImplementation
%
% A more stable version of parametric active contour ("Snake") for safe path planning.
% Key modifications:
%   1) obstacle interior cost = 50 (instead of 1e5),
%   2) cost map = 1 + alphaSoft * exp(-betaSoft*dist), so gradient is gentler,
%   3) bigger internal-energy weights (alpha=5, beta=2),
%   4) smaller dt=1e-4 => slower, more stable updates,
%   5) arc-length reparam every iteration (reparamFreq=1),
%   6) domain clamp + obstacle clamp:
%        if cMap(...) > obsThreshold => snap to boundary
%
% Requirements:
%   - environment.m -> obsMap, startPos, goalPos, xvec,yvec, nx,ny
%   - Then run this code.
%

%% A) Check environment
whosEnv = evalin('base','whos');
envVars = {whosEnv.name};
reqVars = {'obsMap','startPos','goalPos','xvec','yvec','nx','ny'};
for rv = reqVars
    if ~ismember(rv{1}, envVars)
        error(['Variable ',rv{1},' not found in base workspace. ', ...
            'Please run environment.m first.']);
    end
end

%% B) Load
obsMap   = evalin('base','obsMap');   % (ny x nx), 1=obstacle, 0=free
startPos = evalin('base','startPos');
goalPos  = evalin('base','goalPos');
xvec     = evalin('base','xvec');
yvec     = evalin('base','yvec');
nx       = evalin('base','nx');
ny       = evalin('base','ny');

%% C) Build a "soft" cost map
distToObs = bwdist(logical(obsMap));

% inside obstacle => cost = 50
% outside => 1 + alphaSoft * exp(-betaSoft * dist)
% reason: 50 is big, but not infinite -> gradient near boundary still nonzero
cMap = ones(ny,nx)*0;  % pre-alloc
cInside = 50;  % big but finite
for i=1:ny
    for j=1:nx
        if obsMap(i,j)==1
            cMap(i,j) = cInside;
        else
            % outside
            d = distToObs(i,j);
            alphaSoft=5; betaSoft=0.3;
            cMap(i,j) = 1 + alphaSoft*exp(-betaSoft*d);
        end
    end
end

% Precompute gradient for external force
[gradX, gradY] = computeGradient2D(cMap, xvec, yvec);

%% D) Snake parameters
N = 80;  % # points => bigger for smoother curve
gamma = zeros(N,2);
for i=1:N
    frac=(i-1)/(N-1);
    gamma(i,:) = (1-frac)*startPos + frac*goalPos;
end

% PDE iteration controls
dt      = 1e-4;   % smaller => slower but more stable
maxIter = 2000;   % might need more iteration
tol     = 1e-5;   % stricter tolerance

% internal energy
alpha=5.0;  % tension
beta =2.0;  % rigidity

% reparam every iteration
reparamFreq=1;

% Obstacle clamp threshold
% If cMap> clampObs => inside obstacle => snap to boundary
clampObs= cInside - 1;  % e.g. 49

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) First figure: real-time iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Snake with Enhanced Stability','Color','w');
imagesc(xvec,yvec,(cMap>=cInside-1)); % black=obs
colormap([1 1 1; 0 0 0]);
set(gca,'YDir','normal'); axis equal tight; hold on;
plot(startPos(1),startPos(2),'bs','MarkerFaceColor','b','MarkerSize',8);
plot(goalPos(1),goalPos(2),'gs','MarkerFaceColor','g','MarkerSize',8);
title('Snake: stable approach (alpha=5, beta=2, dt=1e-4), obs cost=50');
drawnow;

%% E) Main iteration
for iter=1:maxIter
    gammaOld = gamma;

    % 1) second derivative
    d2 = zeros(N,2);
    for i=2:N-1
        d2(i,:) = gamma(i-1,:) - 2*gamma(i,:) + gamma(i+1,:);
    end

    % 2) fourth derivative
    d4 = zeros(N,2);
    for i=2:N-1
        d4(i,:) = d2(i-1,:) - 2*d2(i,:) + d2(i+1,:);
    end

    % 3) external force
    fext = zeros(N,2);
    for i=1:N
        xNow = gamma(i,1);
        yNow = gamma(i,2);
        [~, gx, gy] = sampleCostClamp(xNow,yNow,cMap,gradX,gradY,xvec,yvec);
        fext(i,:) = -[gx, gy];
    end

    % 4) PDE update
    update = alpha*d2 + beta*d4 + fext;

    for i=2:N-1
        gamma(i,:) = gamma(i,:) + dt*update(i,:);
    end

    % domain clamp
    for i=2:N-1
        gamma(i,1) = max(xvec(1), min(gamma(i,1), xvec(end)));
        gamma(i,2) = max(yvec(1), min(gamma(i,2), yvec(end)));
    end

    % obstacle clamp
    for i=2:N-1
        cVal= sampleCostOnly(gamma(i,1), gamma(i,2), cMap, xvec,yvec);
        if cVal> clampObs
            % snap to obstacle boundary => find local min dist => naive approach
            % for simplicity, let's just shift it a bit outward in direction of negative grad c
            [~,gx,gy]= sampleCostClamp(gamma(i,1),gamma(i,2), cMap, gradX,gradY,xvec,yvec);
            dir = [gx,gy]; 
            nrm= norm(dir);
            if nrm>1e-12
                dir= dir/nrm;
                % move point by small step outward
                step= 0.2; % might need tuning
                gamma(i,:) = gamma(i,:) - step*dir;  
            end
            % also clamp domain
            gamma(i,1) = max(xvec(1), min(gamma(i,1), xvec(end)));
            gamma(i,2) = max(yvec(1), min(gamma(i,2), yvec(end)));
        end
    end

    % fix boundary
    gamma(1,:)   = startPos;
    gamma(N,:)   = goalPos;

    % reparam each iteration
    gamma = reparamArclength(gamma);

    % check difference
    diffMax = max(vecnorm(gamma - gammaOld,2,2));
    if diffMax<tol
        fprintf('Converged at iter=%d, diff=%.4g\n', iter,diffMax);
        break;
    end

    % real-time plot occasionally
    if mod(iter,200)==0 || iter==1
        plot(gamma(:,1), gamma(:,2),'r-','LineWidth',1);
        drawnow;
    end
end

% final path
plot(gamma(:,1), gamma(:,2),'ro-','LineWidth',2,'MarkerFaceColor','r');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Second figure: final path only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Snake Final Path','Color','w');
imagesc(xvec,yvec,(cMap>=cInside-1));
colormap([1 1 1;0 0 0]); 
set(gca,'YDir','normal'); axis equal tight; hold on;
plot(startPos(1),startPos(2),'bs','MarkerFaceColor','b','MarkerSize',8);
plot(goalPos(1),goalPos(2),'gs','MarkerFaceColor','g','MarkerSize',8);
plot(gamma(:,1),gamma(:,2),'ro-','LineWidth',2,'MarkerFaceColor','r');
title('Snake Final Path Only');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER: domain clamp + cost gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cVal,gxVal,gyVal] = sampleCostClamp(x,y,cMap,gradX,gradY,xv,yv)
xC = max(xv(1), min(x, xv(end)));
yC = max(yv(1), min(y, yv(end)));

cVal= interp2(xv,yv,cMap, xC,yC,'linear',50); 
gxVal=interp2(xv,yv,gradX,xC,yC,'linear',0);
gyVal=interp2(xv,yv,gradY,xC,yC,'linear',0);
end

function cVal= sampleCostOnly(x,y,cMap,xv,yv)
xC = max(xv(1), min(x, xv(end)));
yC = max(yv(1), min(y, yv(end)));
cVal= interp2(xv,yv,cMap, xC,yC,'linear',50);
end

function [gradX,gradY] = computeGradient2D(cMap, xvec,yvec)
[ny,nx]=size(cMap);
gradX=zeros(ny,nx);
gradY=zeros(ny,nx);
dx= xvec(2)-xvec(1);
dy= yvec(2)-yvec(1);

for j=2:nx-1
    gradX(:,j)= (cMap(:,j+1)- cMap(:,j-1)) /(2*dx);
end
for i=2:ny-1
    gradY(i,:) = (cMap(i+1,:) - cMap(i-1,:)) /(2*dy);
end
gradX(:,1)= (cMap(:,2)-cMap(:,1))/dx;
gradX(:,nx)=(cMap(:,nx)-cMap(:,nx-1))/dx;
gradY(1,:)=(cMap(2,:)- cMap(1,:))/dy;
gradY(ny,:)=(cMap(ny,:)-cMap(ny-1,:))/dy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER: Reparameterize the snake by arc-length %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gammaOut = reparamArclength(gammaIn)
N=size(gammaIn,1);
segLen=zeros(N-1,1);
for i=1:N-1
    segLen(i)= norm(gammaIn(i+1,:)- gammaIn(i,:));
end
Ltot= sum(segLen);
if Ltot<1e-12
    gammaOut= gammaIn; return;
end
cumLen= [0; cumsum(segLen)];
desired= linspace(0,Ltot,N);
gammaOut=zeros(N,2);
gammaOut(1,:)= gammaIn(1,:);
gammaOut(N,:)= gammaIn(end,:);
idx=1;
for i=2:N-1
    distWanted= desired(i);
    while cumLen(idx+1)< distWanted
        idx= idx+1;
        if idx>=length(cumLen), break; end
    end
    ratio= (distWanted- cumLen(idx)) /(cumLen(idx+1)- cumLen(idx));
    pA= gammaIn(idx,:);
    pB= gammaIn(idx+1,:);
    gammaOut(i,:)= pA + ratio*(pB- pA);
end
end
