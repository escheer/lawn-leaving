% % % % % % 01-06-2017
% % % % % % Elias Scheer
% % % % % %
% % % % % % with tracking code copied from
% % % % % % \MATLAB2016b\toolbox\vision\visiondemos\multiObjectTracking.m
% % % % % %
% % % % objective: count worms as they enter target lawn(s)
% % % % function [tracks, tracksThatLeft, enter_events, exit_events, blobs_in, blobs_out, aversion_ratio] =  tracking_enter_exit(lawn_string, number_of_frames, level, showtracking_flag, singleworm, outer_crop)

clear; clc; close all;
warnid = 'images:initSize:adjustingMag';
warning('off',warnid);

lawn_string = 'nov15_2017_analyzed';
maxwormnumber = 1;
highmag = 1;
outer_crop = 1;
getpixels = 0;
showtracking_flag = 1;

if highmag
    level = 0.93;
else
    level = 0.85;
end
%% Create the video object
close all;

if highmag %javier's setups make multiple movie files per experiment
    [~, pathname, ~] = uigetfile({'*'}); %this is nice because you can see which ones are in progress
    cd(pathname);
    movielist = dir('*.avi');
    movienames = {movielist.name}';
    
    %     Get next file and delete it from the list of movies to be processed
    [ movienames, videoname ] = getnextvideo( movienames );
    v = VideoReader([pathname videoname]);
    number_of_frames = v.NumberOfFrames;
else
    [videoname, pathname, ~] = uigetfile({'*'}); %this is nice because you can see which ones are in progress
    fileID = fopen([videoname 'INPROGRESS'],'w');
    fclose(fileID);
    v = VideoReader([pathname videoname]);
    number_of_frames = v.NumberOfFrames;
end

close all;

%% Pick event horizon and outer boundary
last_frame = imadjust(rgb2gray(im2double(read(v,number_of_frames))));

if getpixels
    pixpermm = getpixpermm( last_frame, 10 );
else
    pixpermm = 112;
end

if outer_crop
    [ outer_boundary_mask, outer_boundary_line, region_rounded ] = get_outer_boundary( last_frame );
    outer_boundary_line = [outer_boundary_line(:,1)-region_rounded(1) outer_boundary_line(:,2)-region_rounded(2)];%shift to the correct coordinates
else
    width = size(last_frame,2);
    height = size(last_frame,1);
    outer_boundary = imrect(ax2,[1 1 width-1 height-1]);
    outer_boundary_line = bbox2points(getPosition(outer_boundary));
    outer_boundary_mask = (outer_boundary.createMask());
    region_rounded = round(getPosition(outer_boundary));
end
%now mask and crop before getting the event_horizon coordinates
masked = last_frame.*outer_boundary_mask;
cropped = imcrop(masked,region_rounded);

%% Generate background and clean it
[orig_background, clean_background, bgthresh] = getbackground( v, region_rounded, outer_boundary_mask, outer_boundary_line, NaN, 0 );

%% Get the event horizon from the background image
[ev_ho_x, ev_ho_y] = geteventhorizon(orig_background,NaN,NaN); %for the fist instance, there is no old ev_ho yet, give it NaN

imshow(orig_background); hold on;
plot(ev_ho_x,ev_ho_y,'g');
pause();
close all;
%% Count worms in the masked area by background subtraction and thresholding
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
    bgsub = (imcomplement(cropped - clean_background));
    thresh1 = im2bw(bgsub,level);
    imshow(thresh1);
    
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
        bgsub = (imcomplement(cropped - clean_background));
        [~, level, ~] = choose_thresh(bgsub,level);
        close all;
        thresh1 = im2bw(bgsub,level);
        imshow(thresh1);
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

%% TRACKING!
hblob = vision.BlobAnalysis;
hblob.Connectivity = 8;
hblob.ExcludeBorderBlobs = true;

if highmag %javier's setup
    hblob = vision.BlobAnalysis;
    hblob.Connectivity = 8;
    hblob.ExcludeBorderBlobs = true;
    hblob.MaximumCount = maxwormnumber;
    hblob.MinimumBlobArea = 300; %was 400 -- keeping this lower may improve detection when the worm background subtraction is problematic
    hblob.MaximumBlobArea = 2500;
    hblob.OutputDataType = 'double';
else
    hblob.MaximumCount = maxwormnumber;
    hblob.MinimumBlobArea = 400; %was 70 on 12 mm
    hblob.MaximumBlobArea = 1200; %was 300 on 12 mm
    % hblob.MaximumBlobArea = 20000; %extremely temporary! for iDU!
end

% fastest way to read in the movie
v2 = VideoReader([pathname videoname]); %have to re-instantiate this object so we can use the faster readFrame method.
v2.CurrentTime = 0; % start at the beginning movie THIS DOESNT EXIST IN 2014a

%for playback

trackPlayer = vision.VideoPlayer('Position', [1043, 32, 816, 593]);

% initialize tracks
[tracks,tracksThatLeft] = initializeTracks_wholevideo();
nextId = 1;
curr_frame = 1;
videoframe = 1;
% showtracking_flag = 0;
skip_loop = 0;

tic;
while hasFrame(v2) || ~isempty(movienames) %MAIN LOOP
    if ~hasFrame(v2) %change the video here
        disp('CHANGE VIDEO!');
        videoframe = 1;
        [ movienames, videoname ] = getnextvideo( movienames );
        v = VideoReader([pathname videoname]);%just for info
        number_of_frames = v.NumberOfFrames;
        [orig_background, clean_background, bgthresh] = getbackground( v, region_rounded, outer_boundary_mask, outer_boundary_line, bgthresh, 1 );
        [ev_ho_x, ev_ho_y] = geteventhorizon(orig_background, ev_ho_x, ev_ho_y);%get the event horizon again
%         figure();
%         imshow(orig_background); hold on;
%         plot(ev_ho_x,ev_ho_y,'g');
%         title(videoname);
        clear v2;
        v2 = VideoReader([pathname videoname]);%for tracking
        v2.CurrentTime = 0;
    end
    disp(['curr frame: ' num2str(curr_frame)]);
    
    frame = imcomplement(im2double(rgb2gray(readFrame(v2))));
    masked = frame.*outer_boundary_mask;
    cropped = imcrop(masked,region_rounded);
    bgsub = (imcomplement(cropped - clean_background));
    thresh1 = imcomplement(im2bw(bgsub,level));
    thresh_cleaned1 = bwareaopen(thresh1,250); %This value was 150
    BWfinal1 = thresh_cleanup( thresh_cleaned1, 29, 3 );
    
    %%% Detect worms: %%%
    [area, centroids, bboxes] = step(hblob,BWfinal1); %image coordinates
    area = double(area);
    bboxes = double(bboxes);
    %%%%%%%%%%%%%%%%%%%%%
    cropworms = cell(size(bboxes,1),1);
    cropworms_orig = cell(size(bboxes,1),1);
    splines = cell(size(bboxes,1),1);
    end1s = NaN(size(bboxes,1),2);
    end2s = NaN(size(bboxes,1),2);
    curvatures = cell(size(bboxes,1),1);
    posture_angles = cell(size(bboxes,1),1);
    worms = cell(size(bboxes,1),1);
    
    for j = 1:size(bboxes,1)
        bbox_crop = [bboxes(j,1)-5 bboxes(j,2)-5 bboxes(j,3)+10 bboxes(j,4)+10];%provide some pixel padding
        wormcrop = imcrop(BWfinal1, bbox_crop);
        wormcrop = imclearborder(wormcrop, 4);
        cropworms{j} = wormcrop;
        wormcrop_orig = imcomplement(imcrop(cropped,bbox_crop)); %would it be better as just cropped instead of bgsub?
        cropworms_orig{j} = wormcrop_orig;
        %only calculate the spline if the centroid is a minimum
        %distance away from the edge of the video
        [~,dist,~] = distance2curve(outer_boundary_line,centroids(j,:));
        if dist > 30 %allowed to get spline at this distance!
            %using Ev's code
            try
                [worm,vWorm,errNum,errMsg] = segWorm_Elias(wormcrop,wormcrop_orig,curr_frame, 0, 1);
                if isempty(errNum) %everything's groovy
                    pt1 = worm.skeleton.pixels(1,:);
                    pt1 = [pt1(2)+bbox_crop(1) pt1(1)+bbox_crop(2)];
                    pt2 = worm.skeleton.pixels(end,:);
                    pt2 = [pt2(2)+bbox_crop(1) pt2(1)+bbox_crop(2)];
                elseif ~isempty(worm) && errNum~= 101 && errNum ~= 102 && errNum ~= 103 && errNum ~= 105 %an error was raised but it's still fine
                    pt1 = worm.skeleton.pixels(1,:);
                    pt1 = [pt1(2)+bbox_crop(1) pt1(1)+bbox_crop(2)];
                    pt2 = worm.skeleton.pixels(end,:);
                    pt2 = [pt2(2)+bbox_crop(1) pt2(1)+bbox_crop(2)];
                elseif isempty(worm) && ~isempty(vWorm.skeleton.pixels) %only save the endpoints, these can be helpful
                    worm = Inf;
                    pt1 = vWorm.skeleton.pixels(1,:);
                    pt1 = [pt1(2)+bbox_crop(1) pt1(1)+bbox_crop(2)];
                    pt2 = vWorm.skeleton.pixels(end,:);
                    pt2 = [pt2(2)+bbox_crop(1) pt2(1)+bbox_crop(2)];
                else %the worm could not be segmented at all
                    worm = Inf;
                    pt1 = [NaN NaN];
                    pt2 = [NaN NaN];
                end
            catch
                worm = Inf;
                pt1 = [NaN NaN];
                pt2 = [NaN NaN];
            end
            skip_loop = 0;
        else
            disp('TOO CLOSE TO BOUNDARY');
            skip_loop = 1;
            break;
        end
        worms{j} = worm;
        if isstruct(worm)
            [spline, curvature, posture_angle] = getWormSpline2( worm, bbox_crop(1), bbox_crop(2) );
        else
            spline = NaN;
            curvature = NaN;
            posture_angle = NaN;
        end
        splines{j} = spline;
        curvatures{j} = curvature;
        posture_angles{j} = posture_angle;
        end1s(j,:) = pt1;
        end2s(j,:) = pt2;
    end
    
    if skip_loop
        if ~isempty(tracks) %update the consecutive invisible counts of tracks in this case
            disp('TRACKS EXIST AND WILL END');
            unassignedTracks = 1:length(tracks); %indices
            [tracks, tracksThatLeft] = updateUnassignedTracks_wholevideo(tracks,tracksThatLeft,unassignedTracks);
        end
        curr_frame = curr_frame + 1;
        videoframe = videoframe +1;
        continue;
    end
    
    tracks = predictNewLocationsOfTracks_wholevideo(tracks);
    
    [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment_wholevideo(tracks,centroids,highmag);
    
    if length(unassignedTracks) >= 1
        warning('WE LOST THE TRACK (NOT OUT OF BOUNDS)');
        disp(curr_frame);
    end
    
    tracks = updateAssignedTracks_wholevideo(tracks, assignments, ev_ho_x, ev_ho_y, centroids, bboxes, cropworms, cropworms_orig, splines, end1s, end2s, worms, curvatures, posture_angles, curr_frame, videoframe, videoname);
    
    [tracks, tracksThatLeft] = updateUnassignedTracks_wholevideo(tracks,tracksThatLeft,unassignedTracks);
    
    [tracks, nextId] = createNewTracks_wholevideo(tracks,nextId,ev_ho_x,ev_ho_y,unassignedDetections,centroids,bboxes,cropworms,cropworms_orig,splines,end1s, end2s, worms,curvatures, posture_angles, curr_frame,videoframe,videoname);
    
    if showtracking_flag
        displayTrackingResults_wholevideo_lawn(trackPlayer, tracks, BWfinal1, ev_ho_x, ev_ho_y);
    end
    curr_frame = curr_frame + 1;
    videoframe = videoframe +1;
end
toc;

%%

%DO POST-PROCESSING ON TRACKS TO DERIVE SPEED, ANGULAR VELOCITY --> ROAMING / DWELLING, LAWN LEAVING EVENTS
allTracks = [tracksThatLeft,tracks];
allTracks = tracks_postprocessing( allTracks, region_rounded(1), region_rounded(2), pixpermm );
%%
% fields2remove = {'cropworm','cropworm_orig','worm','enter_events','exit_events'};
fields2remove = {'cropworm','cropworm_orig','worm'}; %<--- USE THIS ONE
allTracks_slim = rmfield(allTracks,fields2remove);

% IN THE FUTURE
save([videoname '_' lawn_string '_SLIM' '.mat'],'allTracks_slim');
save([videoname '_' lawn_string '.mat'],'allTracks','-v7.3');

% end


