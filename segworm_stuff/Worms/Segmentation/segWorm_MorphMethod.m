function [skel, pt1, pt2, omega] = segWorm_MorphMethod(img, x_offset, y_offset)
skel = NaN; %defaults
pt1 = [NaN NaN];
pt2 = [NaN NaN];
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

% No worm found.
if isempty(wormPixels)
    error('No worm was found.');
end

% Find a point on the contour.
[y, x] = ind2sub(size(img), min(wormPixels));

% Trace the contour clockwise.
contour = bwClockTrace(img, [x y], true);

% The worm is roughly divided into 24 segments of musculature (i.e., hinges
% that represent degrees of freedom) on each side. Therefore, 48 segments
% around a 2-D contour.
% Note: "In C. elegans the 95 rhomboid-shaped body wall muscle cells are
% arranged as staggered pairs in four longitudinal bundles located in four
% quadrants. Three of these bundles (DL, DR, VR) contain 24 cells each,
% whereas VL bundle contains 23 cells." - www.wormatlas.org
sWormSegs = 24;
cWormSegs = 2 * sWormSegs;

cleanContour = cleanWorm(contour, size(contour, 1) / cWormSegs);

%Smooth the contour with a moving mean filter. %ELIAS ADDITION
windowsize = 5;
smthContour = movmean(cleanContour,windowsize,1);

% Compute the contour's local high/low-frequency curvature.
% Note: worm body muscles are arranged and innervated as staggered pairs.
% Therefore, 2 segments have one theoretical degree of freedom (i.e. one
% approximation of a hinge). In the head, muscles are innervated
% individually. Therefore, we sample the worm head's curvature at twice the
% frequency of its body.
% Note 2: we ignore Nyquist sampling theorem (sampling at twice the
% frequency) since the worm's cuticle constrains its mobility and practical
% degrees of freedom.
cCCLengths = circComputeChainCodeLengths(smthContour);
wormSegLength = (cCCLengths(1) + cCCLengths(end)) / cWormSegs;
hfAngleEdgeLength = wormSegLength;
hfCAngles = circCurvature(smthContour, hfAngleEdgeLength, cCCLengths);
lfAngleEdgeLength = 2 * hfAngleEdgeLength;
lfCAngles = circCurvature(smthContour, lfAngleEdgeLength, cCCLengths);

% Blur the contour's local high-frequency curvature.
% Note: on a small scale, noise causes contour imperfections that shift an
% angle from its correct location. Therefore, blurring angles by averaging
% them with their neighbors can localize them better.
wormSegSize = size(smthContour, 1) / cWormSegs;
hfAngleEdgeSize = wormSegSize;
hfBlurSize = ceil(hfAngleEdgeSize / 2);
hfBlurWin(1:hfBlurSize) = 1 / hfBlurSize;
mhfCAngles = circConv(hfCAngles, hfBlurWin);

% Compute the contour's local high/low-frequency curvature maxima.
[mhfCMaxP, ~] = maxPeaksCircDist(mhfCAngles, hfAngleEdgeLength, cCCLengths);
% [lfCMaxP, lfCMaxI] = maxPeaksCircDist(lfCAngles, lfAngleEdgeLength, cCCLengths);

% Check if the worm is doing an omega; this is an unsegmentable posture
% Are the head and tail on the outer contour?
mhfHT = mhfCMaxP > 50;
%ELIAS CHANGED THIS THRESHOLD 11/08/2017 to 45
%ELIAS CHANGED THIS THRESHOLD 01/04/2019 to 50
mhfHTSize = sum(mhfHT);
if mhfHTSize < 2
    disp('Omega! Not segmented.');
    omega = true;
    return;
end

%%% Find properties...
BW_thin             = bwmorph(img,'thin',inf);
BW_endpoints        = bwmorph(BW_thin,'endpoints',inf);
% End points
[I_all ,J_all]      = find(BW_endpoints==1);
NumOfEndPoints      = length(I_all);
if NumOfEndPoints == 0       % e.g. for a circular worm object
    return;
end

branchpoints = bwmorph(BW_thin,'branchpoints');
[I_all, J_all] = find(branchpoints==1);
NumBranchPoints = length(I_all);
% Correct Skeleton if more than 2 endpoints were found: (this method
% chooses the line segment which results in less curvature
% discontinuity) METHOD 1 Elias
BW_pruned = BW_thin;
totalNumPixelsBeforePruning = sum(sum(BW_thin));
counter = 0;
while NumBranchPoints>0
    counter = counter+1;
    bpoint = [I_all(1) J_all(1)];
    %remove branchpoint and all of its 4-connected neighbors to reveal
    %line segments emanating from it.
    north = [bpoint(1)-1 bpoint(2)]; south = [bpoint(1)+1 bpoint(2)];
    east = [bpoint(1) bpoint(2)+1]; west = [bpoint(1) bpoint(2)-1];
    pixelstoremove = [bpoint;north;south;east;west];
    BW_pruned(sub2ind(size(img),pixelstoremove(:,1),pixelstoremove(:,2)))=0; %remove this branchpoint to reveal line segments emanating from it.
    segments = bwconncomp(BW_pruned);
    arclengths = zeros(segments.NumObjects,1);
    for k = 1:segments.NumObjects
        arclengths(k) = length(segments.PixelIdxList{k});
    end
    [~,longestSegIdx] = max(arclengths); %get the longest segment
    [longestSegPixelsI,longestSegPixelsJ] = ind2sub(size(img),segments.PixelIdxList{longestSegIdx});
    %check which concatenated spline results in the least curvature
    %discontinuity
    otherSegs = setdiff(1:segments.NumObjects,longestSegIdx);
    %         maxCurv = Inf;
    maxProm = Inf;
    for k = 1:length(otherSegs)
        otherSegPixels = segments.PixelIdxList{otherSegs(k)};
        [otherSegI,otherSegJ] = ind2sub(size(img),otherSegPixels);
        allPixelstoTest = [longestSegPixelsI longestSegPixelsJ; bpoint; otherSegI otherSegJ]; %concatenate
        tmpImg = zeros(size(img));
        tmpImg(sub2ind(size(img),allPixelstoTest(:,1),allPixelstoTest(:,2)))=1; %set these pixels to 1
        tmpImg = bwmorph(tmpImg,'bridge');
        tmpImg = bwmorph(tmpImg,'thin',inf);
        %to trace boundary of this skeleton in order, get endpoints first
        tmpEndpoints = bwmorph(tmpImg,'endpoints');
        endpts = find(tmpEndpoints==1);
        [eI, eJ] = ind2sub(size(img),endpts);
        sortedConcatSegs = bwtraceboundary(tmpImg, [eI(1) eJ(1)], 'NE'); %get these pixels in order
        Vertices = sortedConcatSegs(1:ceil((end+1)/2),:);
        if length(Vertices)< 0.6*totalNumPixelsBeforePruning %if the current spline greatly reduces the length of spline, it's wrong
            continue;
        end
        BPindex = dsearchn(Vertices,bpoint); %the index of the closest point in spline to the branchpoint
        smthVerts = [movmean(Vertices(:,1),3) movmean(Vertices(:,2),3)];
        curv = CurvatureWindow_Elias(smthVerts, hfAngleEdgeLength*2.5);
        [~,locs,~,proms] = findpeaks(abs(curv));
        peakdist = sqrt((locs-BPindex).^2);
        [closestpeakdist,closestpeakidx] = min(peakdist);
        prom=0;
        if closestpeakdist<5
            prom = proms(closestpeakidx);
        end
        if prom<maxProm
            maxProm = prom;
            BW_pruned = tmpImg;
        end
        %             CurvatBP = abs(curv(BPindex));
        %             if CurvatBP < maxCurv
        %                 maxCurv = CurvatBP;
        %                 BW_pruned = tmpImg;
        %             end
        
    end
    %delete the line segment that resulted in more discontinuous
    %curvature, then check for more branchpoints
    branchpoints = bwmorph(BW_pruned,'branchpoints');
    [I_all, J_all] = find(branchpoints==1);
    NumBranchPoints = length(I_all);
    if counter>100
        disp(['Stopping while loop for skeleton correction', char(10)]);
        break;
    end
end
cc = bwconncomp(BW_pruned);
if cc.NumObjects>1 %sometimes this procedure ends up with disconnected skeletons, this is unsegmentable
    return;
end
if sum(sum(BW_pruned))< 0.6*totalNumPixelsBeforePruning %as before, if the segmentation is too short, discard.
    return;
end
BW_thin = BW_pruned;

% from Hums 2016 code in WormAngles\steps\FindAllWormSkeletonsFromTracks.m
Endpoints = bwmorph(BW_thin, 'endpoints');
[EndpontIndcsX, EndpontIndcsY] = find(Endpoints==1);
if EndpontIndcsX
    wrm = bwtraceboundary(BW_thin, [EndpontIndcsX(1) EndpontIndcsY(1)], 'NE');
    %partition into subsegments
    if (wrm)
        wrmtr=wrm(1:ceil((end+1)/2),:);
        smoothwormX = wrmtr(:,2)';
        smoothwormY = wrmtr(:,1)';
        if numel(smoothwormY) > 45 %commands below will cause a crash when skeleton shrinks too much
            skel = [smoothwormX+x_offset; smoothwormY+y_offset]'; %check that you're adding the right one to each.
            pt1 = round(skel(1,:));
            pt2 = round(skel(end,:));
            
            %DEBUG visualization
%             figure(); imshow(img); hold on;
%             plot(smoothwormX,smoothwormY);
%             pause();
%             close();
            
            return;
        else
            disp('the resulting skeleton is too short!');
            return;
        end
    else
        disp('could not trace the worm!');
        return;
    end
else
    disp('the resulting image has no endpoints!');
    return;
end
end
