function [ ev_ho_x, ev_ho_y, lawn_limit_line] = geteventhorizon_rough(th,background,lawn_limit_line)
%GETEVENTHORIZON_ROUGH.M this function takes in the background image (not
%blurred) and extracts the convex hull -- the limits of the lawn boundary,
%not taking into consideration all of the minor variations along the line.

blurbg = imgaussfilt(background,3);

BW = edge(blurbg,'sobel',th,'nothinning');
BW2 = bwareaopen(BW,500); %get rid of schmutz
% BWnobord = imclearborder(BW2, 8);
BWnobord = BW2; %the border clearing operation can sometimes be problematic if there are streaks
if nargin == 2 %this is the first time we are doing this, so select the in-bounds region for event horizon detection
    disp('SELECT LIMITS ON LAWN BOUNDARY...');
    clf;
    ax2 = axes();
    imshow(BWnobord); hold on;
    lawn_limits = imellipse(ax2); wait(lawn_limits); 
    lawn_limit_line = getVertices(lawn_limits);
    lawn_limit_mask = poly2mask(lawn_limit_line(:,1),lawn_limit_line(:,2),size(BWnobord,1),size(BWnobord,2));
    close all;
else
    lawn_limit_mask = poly2mask(lawn_limit_line(:,1),lawn_limit_line(:,2),size(BWnobord,1),size(BWnobord,2));
end
BWdetect = BWnobord.*lawn_limit_mask;

%find thresholded points and take the convex hull
[lawn_y,lawn_x] = find(BWdetect);

K = convhull(lawn_x,lawn_y);

ev_ho_x = lawn_x(K);
ev_ho_y = lawn_y(K);

end