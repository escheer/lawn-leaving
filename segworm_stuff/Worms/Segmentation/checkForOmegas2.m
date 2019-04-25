function omega = checkForOmegas2(img)

%Based on segWorm_Elias.m in the same folder, but rewritten using mostly
%native MATLAB functions.
%01/04/2019

omega = false;

% Find the worm.
cc = bwconncomp(img);
wormPixels = [];
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
    wormPixels = cc.PixelIdxList{maxCCIdx};
end
% % Find the aspect ratio of the worm's bounding box.
% stats = regionprops(cc);
% dim1 = stats.BoundingBox(3); dim2 = stats.BoundingBox(4);
% aspectRatio = max(dim1,dim2)/min(dim1,dim2);

% Find a point on the contour.
[y, x] = ind2sub(size(img), min(wormPixels));

% Trace the contour.
contour = bwtraceboundary(img, [y x], 'NE');

%Smooth the contour with a moving mean filter.
windowsize = 5;
smthContour = movmean(contour,windowsize,1);

%Calculate the curvature of this smoothed contour.
Vertices = smthContour;
Lines=[(1:size(Vertices,1))' (2:size(Vertices,1)+1)']; Lines(end,2)=1;
smthCurvature = abs(LineCurvature2D(Vertices,Lines));%only care about absolute curvature

%Find curvature peaks that are separated by enough distance along the
%contour to be head and tail

curv_thresh = 0.3; % to designate a head or tail point, the curvature should be greater than 

sWormSegs = 24; %this part is from Ev's code
cWormSegs = 2 * sWormSegs;
cCCLengths = circComputeChainCodeLengths(contour);
wormSegLength = (cCCLengths(1) + cCCLengths(end)) / cWormSegs;
MinDist = wormSegLength*6;

% [pks,~] = findpeaks(smthCurvature,'MinPeakHeight',curv_thresh,'MinPeakDistance',MinDist);
[pks, ~] = maxPeaksCircDist(smthCurvature, MinDist, cCCLengths); %this function works for a circular vector

if sum(pks > curv_thresh) < 2
    errNum = 105;
    errMsg = ['The worm contour has less than 2 high-frequency sampled '...
        'convexities (the head and tail). ' ...
        'Therefore, the worm is coiled or obscured and cannot be segmented. NOT USABLE'];
    omega = true;
    return;
end

end
