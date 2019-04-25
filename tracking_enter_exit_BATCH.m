% % % % % % 12-02-17
% % % % % % Elias Scheer
% % % % % %
% % % % % % with some tracking code copied from
% % % % % % \MATLAB2016b\toolbox\vision\visiondemos\multiObjectTracking.m
% % % % % %

function tracking_enter_exit_BATCH(loadconfig)
clc; close all;
warnid = 'images:initSize:adjustingMag';
warning('off',warnid);

lawn_string = 'dec14_2017_batch_analyzed';

getpixels = 0;

level = 0.93;

th = 0.001;

%% Open the directory containing subdirectories of videos
[~, pathname, ~] = uigetfile({'*'}); %this is nice because you can see which ones are in progress
cd(pathname);
files = dir('*');
folders = files([files.isdir]'); %get the subdirectories
folders = folders(~contains({folders.name}','.')); %get rid of the extra folders starting with '.' and '..'

%make a struct for all of the necessary information to process each video:
%movienames,pixpermm
%outer_boundary_line, outer_boundary_mask, region_rounded,
%orig_background, clean_background, level,
%ev_ho_x, ev_ho_y, lawn_limit_mask
outer_directory = pwd;
if ~loadconfig
    batch_struct = struct();
    for i = 1:length(folders) %go into each subdirectory and get all the necessary information for tracking
        curr_dir = folders(i).name;
        cd(curr_dir);
        batch_struct(i).path = pwd;
        movielist = dir('*.avi');
        movienames = {movielist.name}';
        batch_struct(i).movienames = movienames; %save it
        [ movienames, videoname ] = getnextvideo( movienames );
        v = VideoReader(videoname);
        number_of_frames = v.NumberOfFrames;
        
        %Pick event horizon and outer boundary
        last_frame = imadjust(rgb2gray(im2double(read(v,number_of_frames))));
        
        if getpixels
            pixpermm = getpixpermm( last_frame, 10 );
        else
            pixpermm = 112;
        end
        batch_struct(i).pixpermm = pixpermm;
        %get outer boundary mask
        [ outer_boundary_mask, outer_boundary_line, region_rounded ] = get_outer_boundary_rough( last_frame );
        outer_boundary_line = [outer_boundary_line(:,1)-region_rounded(1) outer_boundary_line(:,2)-region_rounded(2)];%shift to the correct coordinates
        batch_struct(i).outer_boundary_mask = outer_boundary_mask;%save it
        batch_struct(i).outer_boundary_line = outer_boundary_line;
        batch_struct(i).region_rounded = region_rounded;
        
        %Generate background and clean it
        
        out_line_adj = [outer_boundary_line(:,1)+region_rounded(1) outer_boundary_line(:,2)+region_rounded(2)];
        [orig_background, clean_background, level, abovethresh] = getbackground3( v, region_rounded, outer_boundary_mask,out_line_adj,level,1, NaN);
        batch_struct(i).orig_background = orig_background;
        batch_struct(i).clean_background = clean_background;
        batch_struct(i).level = level;
        batch_struct(i).abovethresh = abovethresh;
%         batch_struct(i).dwellingworms = dwellingworms;
        %Get the event horizon from the background image
        [ ev_ho_x, ev_ho_y, lawn_limit_mask] = geteventhorizon_rough(th, orig_background);
        batch_struct(i).ev_ho_x = ev_ho_x;
        batch_struct(i).ev_ho_y = ev_ho_y;
        batch_struct(i).lawn_limit_mask = lawn_limit_mask;
        
        close all;
        cd(outer_directory);
        save([videoname(1:end-11) '_' lawn_string '_CONFIG' '.mat'],'batch_struct');
    end
else
    files = dir('*_CONFIG.mat');
    if length(files)==1
        load(files.name);
    else
        error('there are too many CONFIG files!');
    end
end

%% TRACKING!
for vid = 1:length(batch_struct)
    %get info from struct
    cd(batch_struct(vid).path);
    movienames = batch_struct(vid).movienames;
    [ movienames, videoname ] = getnextvideo( movienames );
    v = VideoReader(videoname);
    number_of_frames = v.NumberOfFrames;
    pixpermm = batch_struct(vid).pixpermm;
    
    outer_boundary_mask = batch_struct(vid).outer_boundary_mask;
    outer_boundary_line = batch_struct(vid).outer_boundary_line;
    region_rounded = batch_struct(vid).region_rounded;
    out_line_adj = [outer_boundary_line(:,1)+region_rounded(1) outer_boundary_line(:,2)+region_rounded(2)]; %this is just for background abovethresh detection
    orig_background = batch_struct(vid).orig_background;
    clean_background = batch_struct(vid).clean_background;
    level = batch_struct(vid).level;
    abovethresh = batch_struct(vid).abovethresh;
    ev_ho_x = batch_struct(vid).ev_ho_x;
    ev_ho_y = batch_struct(vid).ev_ho_y;
    lawn_limit_mask = batch_struct(vid).lawn_limit_mask;
    
    %for background debugging
    bg_struct(1).orig_background = orig_background;
    bg_struct(1).clean_background = clean_background;
    bg_struct(1).level = level;
    bg_struct(1).abovethresh = abovethresh;
    
    %set up tracking
    hblob = vision.BlobAnalysis;
    hblob.Connectivity = 8;
    hblob.ExcludeBorderBlobs = true;
    
    hblob = vision.BlobAnalysis;
    hblob.Connectivity = 8;
    hblob.ExcludeBorderBlobs = true;
    hblob.MaximumCount = 3; %allow three blobs, in order to improve tracking
    hblob.MinimumBlobArea = 300; %was 400 -- keeping this lower may improve detection when the worm background subtraction is problematic
    hblob.MaximumBlobArea = 2500;
    hblob.OutputDataType = 'double';
    
    % fastest way to read in the movie
    v2 = VideoReader(videoname); %have to re-instantiate this object so we can use the faster readFrame method.
    v2.CurrentTime = 0; % start at the beginning movie THIS DOESNT EXIST IN 2014a
    
    % initialize tracks
    [tracks,tracksThatLeft] = initializeTracks_wholevideo();
    nextId = 1;
    curr_frame = 1;
    videoframe = 1;
    skip_loop = 0;
    vidindex = 1;
    
    pastframes = 180; %look back 1 minute to find the current grayscale at the head
    gs_lookback = cell(pastframes,1);
    
    tic;
    while hasFrame(v2) || ~isempty(movienames) %MAIN LOOP
        if ~hasFrame(v2) %change the video here
            disp('CHANGE VIDEO!');
            videoframe = 1;
            vidindex = vidindex+1;
            [ movienames, videoname ] = getnextvideo( movienames );
            v = VideoReader(videoname);%just for info
            number_of_frames = v.NumberOfFrames;
            [orig_background, clean_background, level, abovethresh] = getbackground3(v,region_rounded,outer_boundary_mask,out_line_adj,level,0,abovethresh);
            bg_struct(vidindex).orig_background = orig_background;
            bg_struct(vidindex).clean_background = clean_background;
            bg_struct(vidindex).level = level;
            bg_struct(vidindex).abovethresh = abovethresh;
            [ ev_ho_x, ev_ho_y, lawn_limit_mask] = geteventhorizon_rough(th, orig_background,lawn_limit_mask);
            clear v2;
            v2 = VideoReader(videoname);%for tracking
            v2.CurrentTime = 0;
        end
        disp(['curr frame: ' num2str(curr_frame)]);
        
        frame_ui8 = rgb2gray(readFrame(v2)); %uint8
        frame = imcomplement(im2double(frame_ui8));
        %%%%%%%%% mask out abovethresh pixels from the background (added 12/5/17)
        frame_ui8 = regionfill(frame_ui8,abovethresh);
        frame = regionfill(frame,abovethresh);
        %%%%%%%%%
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
        %get a cleaned, worm-removed version of this image for extracting the
        %grayscale value (proxy for bacterial thickness)
        thresh_gs = 0.4;
        cropped_gs = get_clean_grayscale( frame_ui8, region_rounded, outer_boundary_mask, lawn_limit_mask, thresh_gs, BWfinal1, bboxes );
        
        lookback = 0; %boolean -- whether to look back or not
        %this is a first in first out queue
        if curr_frame<=pastframes
            gs_lookback{curr_frame} = cropped_gs;
            pt1_g = NaN;
            pt2_g = NaN;
            lookback = 0;
        else %start looking back
            cropped_gs_detect = gs_lookback{1};
            gs_lookback(1) = [];%remove the first element
            gs_lookback{pastframes} = cropped_gs;%add the new element
            lookback = 1;
        end
        
        cropworms = cell(size(bboxes,1),1);
        cropworms_orig = cell(size(bboxes,1),1);
        splines = cell(size(bboxes,1),1);
        end1s = NaN(size(bboxes,1),2);
        end2s = NaN(size(bboxes,1),2);
        end1s_g = NaN(size(bboxes,1),1);
        end2s_g = NaN(size(bboxes,1),1);
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
                        if lookback
                            pt1_g = cropped_gs_detect(pt1(1),pt1(2)); %get the grayscale value associated with the head and tail
                            pt2_g = cropped_gs_detect(pt2(1),pt2(2));
                        end
                    elseif ~isempty(worm) && errNum~= 101 && errNum ~= 102 && errNum ~= 103 && errNum ~= 105 %an error was raised but it's still fine
                        pt1 = worm.skeleton.pixels(1,:);
                        pt1 = [pt1(2)+bbox_crop(1) pt1(1)+bbox_crop(2)];
                        pt2 = worm.skeleton.pixels(end,:);
                        pt2 = [pt2(2)+bbox_crop(1) pt2(1)+bbox_crop(2)];
                        if lookback
                            pt1_g = cropped_gs_detect(pt1(1),pt1(2)); %get the grayscale value associated with the head and tail
                            pt2_g = cropped_gs_detect(pt2(1),pt2(2));
                        end
                    elseif isempty(worm) && ~isempty(vWorm.skeleton.pixels) %only save the endpoints, these can be helpful
                        worm = Inf;
                        pt1 = vWorm.skeleton.pixels(1,:);
                        pt1 = [pt1(2)+bbox_crop(1) pt1(1)+bbox_crop(2)];
                        pt2 = vWorm.skeleton.pixels(end,:);
                        pt2 = [pt2(2)+bbox_crop(1) pt2(1)+bbox_crop(2)];
                        if lookback
                            pt1_g = cropped_gs_detect(pt1(1),pt1(2)); %get the grayscale value associated with the head and tail
                            pt2_g = cropped_gs_detect(pt2(1),pt2(2));
                        end
                    else %the worm could not be segmented at all
                        worm = Inf;
                        pt1 = [NaN NaN];
                        pt2 = [NaN NaN];
                        pt1_g = NaN;
                        pt2_g = NaN;
                    end
                catch
                    worm = Inf;
                    pt1 = [NaN NaN];
                    pt2 = [NaN NaN];
                    pt1_g = NaN;
                    pt2_g = NaN;
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
            end1s_g(j) = pt1_g;
            end2s_g(j) = pt2_g;
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
        
        [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment_wholevideo(tracks,centroids);
        
        if length(unassignedTracks) >= 1
            warning('WE LOST THE TRACK (NOT OUT OF BOUNDS)');
            disp(curr_frame);
        end
        
        tracks = updateAssignedTracks_wholevideo(tracks, assignments, ev_ho_x, ev_ho_y, centroids, bboxes, cropworms, cropworms_orig, splines, end1s, end2s, end1s_g, end2s_g, worms, curvatures, posture_angles, curr_frame, videoframe, videoname);
        
        [tracks, tracksThatLeft] = updateUnassignedTracks_wholevideo(tracks,tracksThatLeft,unassignedTracks);
        
        [tracks, nextId] = createNewTracks_wholevideo(tracks,nextId,ev_ho_x,ev_ho_y,unassignedDetections,centroids,bboxes,cropworms,cropworms_orig,splines,end1s, end2s, end1s_g, end2s_g, worms,curvatures, posture_angles, curr_frame,videoframe,videoname);
        
        curr_frame = curr_frame + 1;
        videoframe = videoframe +1;
    end
    toc;
    
    %% DO POST-PROCESSING ON TRACKS TO DERIVE SPEED, ANGULAR VELOCITY --> ROAMING / DWELLING, LAWN LEAVING EVENTS
    allTracks = [tracksThatLeft,tracks];
    WorthSavingTracksThatLeft_idx = [allTracks(:).age]'>9; %double check that all tracks are still long enough -- sometimes the last track is too short to be processed
    allTracks = allTracks(WorthSavingTracksThatLeft_idx);
    allTracks = tracks_postprocessing( allTracks, region_rounded(1), region_rounded(2), pixpermm );
    
    %% SAVING!
    fields2remove = {'cropworm','cropworm_orig','worm'}; %<--- USE THIS ONE
    allTracks_slim = rmfield(allTracks,fields2remove);
    
    save([videoname '_' lawn_string '_BACKGROUND.mat'],'bg_struct');
    save([videoname '_' lawn_string '_SLIM' '.mat'],'allTracks_slim');
    save([videoname '_' lawn_string '.mat'],'allTracks','-v7.3');
    
    
end
end
