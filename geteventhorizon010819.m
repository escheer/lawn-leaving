function [ ev_ho_x, ev_ho_y, lawn_limit_line] = geteventhorizon010819(background, outer_boundary_mask, lawn_limit_line)
%GETEVENTHORIZON_ROUGH.M this function takes in the background image (not
%blurred) and extracts the convex hull -- the limits of the lawn boundary,
%not taking into consideration all of the minor variations along the line.

th = 0.003;
blurbg = imgaussfilt(background,2);

BW = edge(blurbg,'sobel',th,'nothinning');
BWnobord = BW.*outer_boundary_mask;
BWclean = bwareaopen(BWnobord,500); %get rid of schmutz

cleanval = 500;
while sum(sum(BWclean))<3000 %if the cleaning operation deletes all pixels, decrement bwareaopen function iteratively
    cleanval = cleanval - 50;
    BWclean = bwareaopen(BWnobord,cleanval);
end

if nargin == 2 %this is the first time we are doing this, so select the in-bounds region for event horizon detection
    se = strel('disk',5);
    BWclose = imclose(BWclean,se);
    BWnoholes = imfill(BWclose,'holes');
    cc = bwconncomp(BWnoholes);
    lawnPixels = [];
    if ~isempty(cc.PixelIdxList)
        maxCCIdx = 0;
        maxCCSize = 0;
        for i = 1:length(cc.PixelIdxList)
            ccSize = length(cc.PixelIdxList{i});
            if ccSize > maxCCSize
                maxCCSize = ccSize;
                maxCCIdx = i;
            end
        end
        lawnPixels = cc.PixelIdxList{maxCCIdx};
    end
    BWdetect = zeros(size(BWnoholes)); %make a mask for the largest connected component
    BWdetect(lawnPixels)=1;
    [contour_x, contour_y] = get_lawn_boundary_from_pixels(BWdetect);
    
    % Fit a circle to this boundary, which will become the lawn-limit line
    Par = CircleFitByPratt([contour_x contour_y]);
    c_x = Par(1); c_y = Par(2); radius = Par(3)+10; %add some to avoid messiness on the boundary
    lawn_limit_line = circlepoints( c_x, c_y, radius );
    figure(); imshow(background); hold on;
    plot(contour_x,contour_y,'LineWidth',2);
    plot(lawn_limit_line(:,1),lawn_limit_line(:,2),'LineWidth',2);
    good = input('LAWN BOUNDARY LOOKS GOOD? (1/0)');
    while ~(good==0 || good == 1)
        good = input('OUTER BOUNDARY LOOKS GOOD? (1/0)');
    end
    close all;
    while ~good %if not good do it the manual way
        disp('MANUAL WAY: SELECT LIMITS ON LAWN BOUNDARY...');
        clf;
        ax2 = axes();
        imshow(BWclean); hold on;
        lawn_limits = imellipse(ax2); wait(lawn_limits);
        lawn_limit_line = getVertices(lawn_limits);
        lawn_limit_mask = poly2mask(lawn_limit_line(:,1),lawn_limit_line(:,2),size(blurbg,1),size(blurbg,2));
        close all;
        BWdetect = BWnobord.*lawn_limit_mask;
        [contour_x, contour_y] = get_lawn_boundary_from_pixels(BWdetect);
        figure(); imshow(background); hold on;
        plot(contour_x,contour_y,'LineWidth',2);
        plot(lawn_limit_line(:,1),lawn_limit_line(:,2),'LineWidth',2);
        good = input('LAWN BOUNDARY LOOKS GOOD? (1/0)');
        while ~(good==0 || good == 1)
            good = input('OUTER BOUNDARY LOOKS GOOD? (1/0)');
        end
    end
    close all;
    ev_ho_x = contour_x;
    ev_ho_y = contour_y;
    
else %upon successive iterations, UI is not required
    lawn_limit_mask = poly2mask(lawn_limit_line(:,1),lawn_limit_line(:,2),size(blurbg,1),size(blurbg,2));
    BWdetect = BWclean.*lawn_limit_mask;
    [ev_ho_x, ev_ho_y] = get_lawn_boundary_from_pixels(BWdetect);
end

end

function [contour_x, contour_y] = get_lawn_boundary_from_pixels(BW) %take in an image with the pixels only
[lawn_y,lawn_x] = find(BW);
K = boundary(lawn_x,lawn_y);
lawnBoundary = [lawn_x(K) lawn_y(K)];

% Smooth lawn contour
smthWindow = 15;
% smthContour = [movmean(lawnBoundary(:,1),smthWindow) movmean(lawnBoundary(:,2),smthWindow)];
smthContour = movmean(lawnBoundary,smthWindow,1,'Endpoints','discard');

% Shift the lawn boundary, and re-smooth, to get rid of gap
lawnBoundary_shift = circshift(lawnBoundary,round(size(lawnBoundary,1)/2),1);
% smthContour_shift = [movmean(lawnBoundary_shift(:,1),smthWindow) movmean(lawnBoundary_shift(:,2),smthWindow)];
smthContour_shift = movmean(lawnBoundary_shift,smthWindow,1,'Endpoints','discard');

missingContour = setdiff(smthContour_shift,smthContour,'rows');
newContour = [smthContour ; zeros(size(missingContour))];
contourEndpoint = smthContour(end,:);
newidx = size(smthContour,1)+1;
while size(missingContour,1)>0
    [~,closeridx] = min(sqrt(sum(bsxfun(@minus, missingContour, contourEndpoint).^2,2))); %find the closest point in the missing contour to add
    newEndpoint = missingContour(closeridx,:);
    newContour(newidx,:) = newEndpoint;
    contourEndpoint = newEndpoint;
    missingContour(closeridx,:)=[];
    newidx = newidx+1;
end
contour_x = newContour(:,1);
contour_y = newContour(:,2);

% % Take the boundary again %THIS METHOD CAN CAUSE SOME IRREGULARITIES.
% bothContour = [smthContour;smthContour_shift];
% K = boundary(bothContour(:,1),bothContour(:,2));
% bothBoundary = [bothContour(K,1) bothContour(K,2)];
% 
% contour_x = bothBoundary(:,1);
% contour_y = bothBoundary(:,2);
end