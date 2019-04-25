function [ outer_boundary_mask, outer_boundary_line, region_rounded ] = get_outer_boundary( frame )
%GET_OUTER_BOUNDARY.M This function estimates the outer boundary line from
%the last frame of the first movie -- this is held as a constant throughout
%the movie

% outer = (edge(frame,'sobel'));
th = 0.1;
outer = (edge(frame,'sobel',th,'nothinning'));
outer = imclearborder(outer,8);

se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
% outer_dil = imdilate(outer,[se90 se0]);

se = strel('disk',5);
outer = imdilate(outer,se);
outer = imdilate(outer,se);
outer = imdilate(outer,se);
outer_dil = imdilate(outer,[se90 se0]);

outer_dil_clean = bwareaopen(outer_dil,10000);
outer_fill = imfill(outer_dil_clean,'holes');

outer_erode = imerode(outer_fill,se);
outer_erode = imerode(outer_erode,se);
outer_erode = imerode(outer_erode,se);
outer_erode = imerode(outer_erode,[se0 se90]);

% outer_clear = imclearborder(outer_erode);
outer_clean = bwareaopen(outer_erode,10000);
pad_outer_clean = padarray(outer_clean,[5 5]);
CC = bwconncomp(pad_outer_clean);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
if ~isempty(idx)
    BWbig = zeros(size(pad_outer_clean));
    BWbig(CC.PixelIdxList{idx})=1;
    
    boundaries = bwboundaries(BWbig);
    
    if size(boundaries,1)==1
        firstBoundary = boundaries{1};
        x = firstBoundary(:, 2)-5; %to correct for the padarray
        y = firstBoundary(:, 1)-5;
        choose_manually = 0;
    else
        choose_manually = 1;
    end
else
    choose_manually = 1;
end

if ~choose_manually %can detect it automatically!
    polynomialOrder = 3;
    windowWidth = 751;        
    smoothX = sgolayfilt(x, polynomialOrder, windowWidth);
    smoothY = sgolayfilt(y, polynomialOrder, windowWidth);
    outer_boundary_line = [smoothX smoothY];
    
    figure(2); ax2 = axes();
    imshow(frame); hold on;
    plot(smoothX, smoothY,'LineWidth',5);
    Par = CircleFitByPratt(outer_boundary_line);
    c_x = Par(1); c_y = Par(2); radius = Par(3)-45; %subtract some to avoid messiness on the boundary
    outer_boundary_line = circlepoints( c_x, c_y, radius );
    plot(outer_boundary_line(:,1),outer_boundary_line(:,2),'LineWidth',5);
    good = input('LOOKS GOOD?');
    if ~good %if not, choose manually
        clf;
        ax2 = axes();
        imshow(frame); hold on;
        outer_boundary = imellipse(ax2); wait(outer_boundary);
        outer_boundary_line = getVertices(outer_boundary);
        plot(outer_boundary_line(:,1),outer_boundary_line(:,2),'LineWidth',5);
        Par = CircleFitByPratt(outer_boundary_line);
        c_x = Par(1); c_y = Par(2); radius = Par(3)-45; %subtract some to avoid messiness on the boundary
        outer_boundary_line = circlepoints( c_x, c_y, radius );
        plot(outer_boundary_line(:,1),outer_boundary_line(:,2),'LineWidth',5);
        pause();
    end
else %choose manually
    figure(2); ax2 = axes();
    imshow(frame); hold on;
    outer_boundary = imellipse(ax2); wait(outer_boundary);
    outer_boundary_line = getVertices(outer_boundary);
    plot(outer_boundary_line(:,1),outer_boundary_line(:,2),'LineWidth',5);
    Par = CircleFitByPratt(outer_boundary_line);
    c_x = Par(1); c_y = Par(2); radius = Par(3)-45; %subtract some to avoid messiness on the boundary
    outer_boundary_line = circlepoints( c_x, c_y, radius );
    plot(outer_boundary_line(:,1),outer_boundary_line(:,2),'LineWidth',5);
end
% Par = CircleFitByPratt(outer_boundary_line);
% c_x = Par(1); c_y = Par(2); radius = Par(3)-45; %subtract some to avoid messiness on the boundary
% outer_boundary_line = circlepoints( c_x, c_y, radius );
% plot(outer_boundary_line(:,1),outer_boundary_line(:,2),'LineWidth',5);
% pause();
outer_boundary_mask = poly2mask(outer_boundary_line(:,1),outer_boundary_line(:,2),size(frame,1),size(frame,2));
RP = regionprops(outer_boundary_mask);
region_rounded = round(RP.BoundingBox);

close all;
end