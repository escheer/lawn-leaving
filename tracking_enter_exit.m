% % % 01-06-2017
% % % Elias Scheer
% % % inspired by code from Javier Marquina-Solis, which was modified from Frank Tejera to count worms.
% % % with tracking code copied from
% % % \MATLAB2016b\toolbox\vision\visiondemos\multiObjectTracking.m
% % %
% objective: count worms as they enter target lawn(s)
function [tracks, tracksThatLeft, enter_events, exit_events, blobs_in, blobs_out, aversion_ratio] =  tracking_enter_exit(lawn_string, last_frame_idx, level, showtracking_flag, singleworm, outer_crop)

%% Create the video object
close all;
[filename, pathname, ~] = uigetfile({'*'}); %this is nice because you can see which ones are in progress
fileID = fopen([filename 'INPROGRESS'],'w');
fclose(fileID);
v = VideoReader([pathname filename]);
number_of_frames = v.NumberOfFrames;

uncropped_vid = [filename(1:end-12) '.avi'];
tmpVid = VideoReader([pathname uncropped_vid]);
f1 = im2double(read(tmpVid,1));
imshow(f1);
maxwormnumber = input('Tell me the number of worms in the video!\n');
close all;
% clear tmpVid;
% clear f1;

%% Pick event horizon and outer boundary
h1 = figure(1); ax1 = axes();
last_frame = im2double(read(v,number_of_frames));
% last_frame = imadjust(rgb2gray(im2double(read(v,number_of_frames))));
imshow(imadjust(rgb2gray(last_frame)));
% imshow(last_frame);
if outer_crop
    outer_boundary = imellipse(ax1); wait(outer_boundary);
    outer_boundary_mask = (outer_boundary.createMask());
    outer_region = getPosition(outer_boundary);
    region_rounded = round(outer_region);
else
    width = size(last_frame,2);
    height = size(last_frame,1);
    outer_boundary = imrect(ax1,[1 1 width-1 height-1]);
    outer_boundary_mask = (outer_boundary.createMask());
    outer_region = getPosition(outer_boundary);
    region_rounded = round(outer_region);
end
close all;
%now mask and crop before getting the event_horizon coordinates
h2 = figure(2); ax2 = axes();
masked = last_frame.*outer_boundary_mask;
cropped = imcrop(masked,region_rounded);

%% Generate background
%frames used to generate the background
bg_frames = linspace(floor(number_of_frames/2),number_of_frames,60);
nf = length(bg_frames);
% Fill the background matrix with zeros
stack_back = zeros(region_rounded(4)+1, region_rounded(3)+1, nf);
% Stack all the frames together
for i = 1:nf
    frame = imcomplement(im2double(rgb2gray(read(v,bg_frames(i))))); % read in frame ITS DOUBLE
    masked = frame.*outer_boundary_mask;
    cropped = imcrop(masked,region_rounded);
    stack_back(:,:,i) = cropped;
end

% Mean of the stack to obtain the background
background = (mean(stack_back,3));
% figure, imshow(background)
% pause(); 
close all;


%% get the event horizon from the background image
backgthresh = multithresh(background,2);
imshow(im2bw(background,backgthresh(2)))
bggood = input('Is thresholded background good? (1/0) \n');
if bggood
    bglevel = backgthresh(2);
else
    imshow(im2bw(background,backgthresh(1)))
    lastchance = input('How about now? (1/0) \n');
    if lastchance
        bglevel = backgthresh(1);
    else
        error('problem getting background threshold!')
    end
end
threshbg = im2bw(background,bglevel);
threshbg = bwareaopen(threshbg,1000);
[B,~,N] = bwboundaries(threshbg,'noholes');
if N > 1
    imshow(threshbg); hold on;
    pause();
    error('there can only be one lawn!');
end
boundary = B{1};
ev_ho_x = boundary(:,2);
ev_ho_y = boundary(:,1);
imshow(background); hold on;
scatter(ev_ho_x,ev_ho_y,'g+');
pause();

% cc = bwconncomp(threshbg);
% if cc.NumObjects > 1
%     imshow(threshbg); hold on;
%     pause();
%     error('there can only be one lawn!');
% end
% cc_obj = cc.PixelIdxList{1}
% [cc_i, cc_j] = ind2sub(cc.ImageSize,cc_obj);
% imshow(threshbg); hold on;
% scatter(cc_j,cc_i,'r+')


% imshow(imadjust(rgb2gray(cropped)));
% event_horizon = imellipse(ax2); wait(event_horizon);
% pos_ev_ho = getPosition(event_horizon);
% %get points along this ellipse for checking whether tracks have passed
% %through it
% xmin = pos_ev_ho(1); ymin = pos_ev_ho(2); width = pos_ev_ho(3); height = pos_ev_ho(4);
% xCenter = (xmin+xmin+width)/2; yCenter = (ymin+ymin+height)/2;
% halfwidth = width/2; halfheight = height/2;
% ev_ho_ellipse = ellipse(halfwidth,halfheight,0,xCenter,yCenter); %make sure this is in the path
% ev_ho_x = ev_ho_ellipse.XData; ev_ho_y = ev_ho_ellipse.YData;

close all;




%% Count worms in the masked area by background subtraction and thresholding
hblob = vision.BlobAnalysis;
hblob.Connectivity = 8;
hblob.ExcludeBorderBlobs = true;

if singleworm %javier's setup
    hblob.MaximumCount = 3;
    hblob.MinimumBlobArea = 400;
    hblob.MaximumBlobArea = 2000;
else
    hblob.MaximumCount = 100;
    hblob.MinimumBlobArea = 80;
    hblob.MaximumBlobArea = 200;
end


% initialize tracks
[tracks,tracksThatLeft] = initializeTracks_wholevideo();
nextId = 1;

% fastest way to read in the movie
v2 = VideoReader([pathname filename]); %have to re-instantiate this object so we can use the faster readFrame method.
v2.CurrentTime = 0; % start at the beginning movie THIS DOESNT EXIST IN
% 2014a

%for playback
if showtracking_flag
    blobPlayer = vision.VideoPlayer('Position', [120, 160, 850, 850]);
    trackPlayer = vision.VideoPlayer('Position', [970, 160, 850, 850]);
end



if nargin<3 %if no level is given, go through a bunch of random frames and average the otsu's level from it
    levels = zeros(100,1);
    frames_to_check = randi(number_of_frames,100,1);
    for i = length(frames_to_check)
        frame = imcomplement(im2double(rgb2gray((read(v,frames_to_check(i))))));
        masked = frame.*outer_boundary_mask;
        cropped = imcrop(masked,region_rounded);
        bgsub = (imcomplement(cropped - background));
        levels(i) = graythresh(bgsub); %otsu's method for getting a threshold level for binarization
    end
    level = mean(levels)+0.1; % I find Otsu's method yields systematically too low thresholds
end

frames_to_check = randi([fix(number_of_frames/2) number_of_frames],100,1);
goodcounter = 0;
badcounter = 0;
for j = 1:length(frames_to_check)
    if goodcounter > 3 || badcounter > 0
        break;
    end
    frame = imcomplement(im2double(rgb2gray((read(v,frames_to_check(j))))));
    masked = frame.*outer_boundary_mask;
    cropped = imcrop(masked,region_rounded);
    bgsub = (imcomplement(cropped - background));
    close all;
    thresh = im2bw(bgsub,level);
    imshow(thresh);
    
    
    ButtonName = questdlg('Is threshold OK?', ...
        'Check Threshold', ...
        'GOOD', 'BAD', 'NO WORMS', 'NO WORMS');
    switch ButtonName,
        case 'BAD',
            badcounter = badcounter+1;
        case 'GOOD',
            goodcounter = goodcounter+1;
        case 'NO WORMS',
            continue;
    end
end
close all;
%if there is something wrong, choose the new threshold
if badcounter > 0
    levels = NaN;
    goodcounter = 0;
    for j = 1:length(frames_to_check)
        if goodcounter > 3
            break;
        end
        frame = imcomplement(im2double(rgb2gray((read(v,frames_to_check(j))))));
        masked = frame.*outer_boundary_mask;
        cropped = imcrop(masked,region_rounded);
        bgsub = (imcomplement(cropped - background));
        [~, level, ~] = choose_thresh(bgsub,level);
        close all;
        thresh = im2bw(bgsub,level);
        imshow(thresh);
        button = questdlg('is threshold good (with worms in it)?');
        if strcmp(button,'No')
            continue;
        else
            levels = [levels level];
            goodcounter = goodcounter+1;
        end
    end
    close all;
    level = nanmean(levels);
end


tic;
loop_counter = 1;

while hasFrame(v2) %main loop

    frame = imcomplement(im2double(rgb2gray(readFrame(v2))));
    masked = frame.*outer_boundary_mask;
    cropped = imcrop(masked,region_rounded);
    bgsub = (imcomplement(cropped - background));
    thresh = imcomplement(im2bw(bgsub,level));
    thresh_cleaned = bwareaopen(thresh,50);
    [~, centroids, bboxes] = step(hblob,thresh_cleaned); %image coordinates
    
    tracks = predictNewLocationsOfTracks_wholevideo(tracks);
    
    [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment_wholevideo(tracks,centroids,singleworm);
    
    tracks = updateAssignedTracks_wholevideo(tracks,assignments,centroids,bboxes,loop_counter);
    
    [tracks, tracksThatLeft] = updateUnassignedTracks_wholevideo(tracks,tracksThatLeft, unassignedTracks,loop_counter);
      
    [tracks, nextId] = createNewTracks_wholevideo(tracks,nextId, unassignedDetections,centroids,bboxes,loop_counter);
    if showtracking_flag
        displayTrackingResults_wholevideo_lawn(blobPlayer, bboxes, trackPlayer, tracks, thresh_cleaned, ev_ho_x, ev_ho_y)
    end
    loop_counter = loop_counter + 1;
    
end
toc;
% [enter_events, exit_events] = countCrossingTracks_after( tracks,tracksThatLeft, ev_ho_x, ev_ho_y );
[blobs_in, blobs_out] = countBlobsInOut(tracks,tracksThatLeft,ev_ho_x, ev_ho_y, last_frame_idx );
aversion_ratio = blobs_out./maxwormnumber;
[enter_events, exit_events] = countCrossingTracks_after( tracks,tracksThatLeft, ev_ho_x, ev_ho_y);

save([filename '_' lawn_string '.mat'],'tracks','tracksThatLeft','enter_events','exit_events','blobs_in','blobs_out','aversion_ratio','ev_ho_x','ev_ho_y');

end


