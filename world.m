%% 1) 도메인 및 그리드 설정
clear; clc; close all;

nx = 300;   % x방향 그리드 갯수
ny = 300;   % y방향 그리드 갯수
xrange = [0, 10];  % x 범위
yrange = [0, 10];  % y 범위

xvec = linspace(xrange(1), xrange(2), nx);
yvec = linspace(yrange(1), yrange(2), ny);
[XX, YY] = meshgrid(xvec, yvec); 
% XX(i,j) = x좌표, YY(i,j) = y좌표

%% 2) 장애물(원형 여러 개) 정의
% [cx, cy, r]
obstacles = [
    2, 4, 2.5;
    8, 2, 2.2;
    7, 7, 1.7;
    2, 9, 2.0
];

% 장애물 지도: obsMap(i,j)=1 -> 장애물, 0 -> 통과 가능
obsMap = zeros(ny, nx);

for k = 1:size(obstacles,1)
    cx = obstacles(k,1);
    cy = obstacles(k,2);
    r  = obstacles(k,3);

    distSq = (XX - cx).^2 + (YY - cy).^2;
    inside = (distSq <= r^2);
    obsMap(inside) = 1;
end

%% 3) 시작점 / 목표점 설정
startPos = [0, 0];
goalPos  = [10, 10];

%% 4) 시각화
figure('Color','w');

% 수동 컬러맵 설정: 
%   0 -> [1 1 1] (흰색 배경)
%   1 -> [1 0 0] (빨간색 장애물)
myColorMap = [1 1 1;    % index=0인 값의 색 (배경)
              0 0 0];   % index=1인 값의 색 (장애물)
colormap(myColorMap);

% obsMap의 값 범위를 0~1로 맞추기 위해 caxis([0 1]) 설정
imagesc(xvec, yvec, obsMap);
caxis([0 1]);
set(gca, 'YDir','normal');
axis equal; axis tight; hold on;
title('Obstacle Map (Obstacle=Red, Free=White)');

% 시작점(파랑), 목표점(초록) 등 눈에 띄는 컬러로 변경
plot(startPos(1), startPos(2), 'bo','MarkerFaceColor','b','MarkerSize',8);
plot(goalPos(1),  goalPos(2),  'go','MarkerFaceColor','r','MarkerSize',8);

% 텍스트 오프셋을 두어 겹치지 않게 하고, 배경색(white), 여백(margin) 설정
textOffset = 0.3;
text(startPos(1)+textOffset, startPos(2)+textOffset, 'Start',...
    'Color','blue','FontWeight','bold',...
    'BackgroundColor','w','Margin',2);
text(goalPos(1)-1.5, goalPos(2)-textOffset, 'Goal',...
    'Color','red','FontWeight','bold',...
    'BackgroundColor','w','Margin',2);

xlim(xrange); ylim(yrange);
