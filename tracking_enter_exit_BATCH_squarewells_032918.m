function tracking_enter_exit_BATCH_012418(loadconfig, getpixels)
clc; close all;
warnid = 'images:initSize:adjustingMag';
warning('off',warnid);

timenow = datestr(now,'mm_dd_yy');
lawn_string = ['012418_BATCH_' timenow];

maxwormnumber = 3;

% getpixels = 0;

level = 0.90;

th = 0.001; %for event horizon detection
%% Open the directory containing subdirectories of videos
[filename, pathname, ~] = uigetfile({'*'}); %this is nice because you can see which ones are in progress
cd(pathname);
files = dir('*');
folders = files([files.isdir]'); %get the subdirectories
folders = folders(~contains({folders.name}','.')); %get rid of the extra folders starting with '.' and '..'

%make a struct for all of the necessary information to process each video

outer_directory = pwd;
if ~loadconfig
    batch_struct = struct();
    %BATCH_STRUCT HAS THE FOLLOWING FIELDS:
    % PATH
    % MOVIENAMES
    % PIXPERMM
    % CROP_POS
    % OUTER_BOUNDARY_MASK
    % OUTER_BOUNDARY_LINE
    % OUTER_BOUNDARY_LINE_CRP_REL
    % REGION_ROUNDED
    % ORIG_BACKGROUND
    % CLEAN_BACKGROUND
    % LEVEL
    % ABOVETHRESH
    % FULLYBLURRED_BACKGROUND
    % LAWN_LIMIT_MASK_WF
    % LANW_LIMIT_LINE
    % EV_HO
    % EV_HO_CRP_REL
    
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
            disp(pixpermm);
        else
            pixpermm = 112;
        end
        batch_struct(i).pixpermm = pixpermm;
        
        % select region for correcting plate bumps
        disp('SELECT REGION OUTSIDE OF THE PLATE FOR FRAME REGISTRATION.');
        x = figure(1);
        ax = axes(); hparent = imshow(last_frame);
        h = imrect(ax,[10 10 100 100]); wait(h);
        crop_pos = getPosition(h);
        close(x);
        batch_struct(i).crop_pos = crop_pos;
        
        %get outer boundary mask
        ax = axes();
        imshow(last_frame); hold on;
        outer_boundary = imrect(ax); wait(outer_boundary);
        outer_boundary_line = round(bbox2points(getPosition(outer_boundary)));
        outer_boundary_mask = poly2mask(outer_boundary_line(:,1),outer_boundary_line(:,2),size(last_frame,1),size(last_frame,2));
        RP = regionprops(outer_boundary_mask);
        region_rounded = round(RP.BoundingBox);
        outer_boundary_line_crp_rel = [outer_boundary_line(:,1)-region_rounded(1) outer_boundary_line(:,2)-region_rounded(2)];%shift to the correct coordinates
        close();
        
        batch_struct(i).outer_boundary_mask = outer_boundary_mask;%save it
        batch_struct(i).outer_boundary_line = outer_boundary_line;
        batch_struct(i).outer_boundary_line_crp_rel = outer_boundary_line_crp_rel;
        batch_struct(i).region_rounded = region_rounded;
        
        %Generate background and clean it
        [orig_background, clean_background, filledbg, level, abovethresh] = getbackground3( v, outer_boundary_mask, outer_boundary_line, level,1, NaN, NaN);
        batch_struct(i).orig_background = orig_background;
        batch_struct(i).clean_background = clean_background;
        batch_struct(i).level = level;
        batch_struct(i).abovethresh = abovethresh;
        
        %Get the event horizon from the background image
        [ ev_ho_x, ev_ho_y, lawn_limit_line] = geteventhorizon_rough(th,orig_background);
        ev_ho_crp_rel = [ev_ho_x-region_rounded(1) ev_ho_y-region_rounded(2)]; %event horizon in cropped coordinates
        lawn_limit_mask = poly2mask(lawn_limit_line(:,1),lawn_limit_line(:,2),size(last_frame,1),size(last_frame,2));
        fullyblurred_background = imgaussfilt(regionfill(filledbg,lawn_limit_mask),4);
        batch_struct(i).fullyblurred_background = fullyblurred_background;
        batch_struct(i).lawn_limit_mask_wf = lawn_limit_mask;  % in whole frame coordinates
        batch_struct(i).lawn_limit_line = lawn_limit_line;     % in whole frame coordinates
        batch_struct(i).ev_ho = [ev_ho_x ev_ho_y];             % in whole frame coordinates
        batch_struct(i).ev_ho_crp_rel = ev_ho_crp_rel;         % in cropped frame coordinates
        
        close all;
        cd(outer_directory);
        save([videoname(1:end-11) '_' lawn_string '_CONFIG' '.mat'],'batch_struct');
    end
else
%     files = dir('*_CONFIG.mat');
%     if length(files)==1
%         load(files.name);
%     else
%         error('there are too many CONFIG files!');
%     end
        load(filename);
end




%% TRACKING!
%BG_STRUCT HAS THE FOLLOWING FIELDS:
% VIDEONAME
% VIDEOFRAME
% CURR_FRAME
% REGION_ROUNDED
% OUTER_BOUNDARY_MASK
% OUTER_BOUNDARY_LINE
% ORIG_BACKGROUND
% CLEAN_BACKGROUND
% LEVEL
% ABOVETHRESH
% FULLYBLURRED_BACKGROUND
% LAWN_LIMIT_MASK_WF
% LAWN_LIMIT_LINE
% EV_HO
% EV_HO_CRP_REL

for vid = 1:length(batch_struct)
    %go to the correct directory
    cd(batch_struct(vid).path);
    %%% GET OUT ALL DATA FROM BATCH STRUCT
    movienames = batch_struct(vid).movienames;
    pixpermm = batch_struct(vid).pixpermm;
    crop_pos = batch_struct(vid).crop_pos;
    outer_boundary_mask = batch_struct(vid).outer_boundary_mask;
    outer_boundary_line = batch_struct(vid).outer_boundary_line;
    outer_boundary_line_crp_rel = batch_struct(vid).outer_boundary_line_crp_rel;
    region_rounded = batch_struct(vid).region_rounded;
    orig_background = batch_struct(vid).orig_background;
    clean_background = batch_struct(vid).clean_background;
    level = batch_struct(vid).level;
    abovethresh = batch_struct(vid).abovethresh;
    fullyblurred_background = batch_struct(vid).fullyblurred_background;
    lawn_limit_mask_wf = batch_struct(vid).lawn_limit_mask_wf;
    lawn_limit_line = batch_struct(vid).lawn_limit_line;
    ev_ho = batch_struct(vid).ev_ho;
    ev_ho_x = ev_ho(:,1); ev_ho_y = ev_ho(:,2); %ADDED 12/25 THIS IS NECESSARY AS THE LOCAL VARIABLE
    ev_ho_crp_rel = batch_struct(vid).ev_ho_crp_rel;
    %%%
    [ movienames, videoname ] = getnextvideo( movienames );
    v = VideoReader(videoname);
    videoframe = 1;
    curr_frame = 1;
    %%% ADD ALL RELEVANT INFORMATION TO THE BG_STRUCT
    bg_struct = struct(); %make it fresh for every new set of videos
    bg_struct(1).videoname = videoname;
    bg_struct(1).videoframe = videoframe;
    bg_struct(1).curr_frame = curr_frame;
    bg_struct(1).region_rounded = region_rounded;
    bg_struct(1).outer_boundary_mask = outer_boundary_mask;
    bg_struct(1).outer_boundary_line = outer_boundary_line;
    bg_struct(1).orig_background = orig_background;
    bg_struct(1).clean_background = clean_background;
    bg_struct(1).level = level;
    bg_struct(1).abovethresh = abovethresh;
    bg_struct(1).fullyblurred_background = fullyblurred_background;
    bg_struct(1).lawn_limit_mask_wf = lawn_limit_mask_wf;
    bg_struct(1).lawn_limit_line = lawn_limit_line;
    bg_struct(1).ev_ho = ev_ho;
    bg_struct(1).ev_ho_crp_rel = ev_ho_crp_rel;
    %%%
    
    %%% SET UP BLOB ANALYSIS
    hblob = vision.BlobAnalysis;
    hblob.Connectivity = 8;
    hblob.ExcludeBorderBlobs = true;
    
    hblob = vision.BlobAnalysis;
    hblob.Connectivity = 8;
    hblob.ExcludeBorderBlobs = true;
    hblob.MaximumCount = maxwormnumber;
    hblob.MinimumBlobArea = 300; %was 400 -- keeping this lower may improve detection when the worm background subtraction is problematic
    hblob.MaximumBlobArea = 2500;
    hblob.OutputDataType = 'double';
    
    % fastest way to read in the movie
    v2 = VideoReader(videoname); %have to re-instantiate this object so we can use the faster readFrame method.
    v2.CurrentTime = 0; % start at the beginning movie THIS DOESNT EXIST IN 2014a
    
    %%% INITIALIZE TRACKS AND VALUES NEEDED FOR TRACKING
    [tracks,tracksThatLeft] = initializeTracks_wholevideo();
    nextId = 1;
    skip_loop = 0;
    bgvidindex = 1;
    pastframes = 360; %look back 2 minute to find the current grayscale at the head
    gs_lookback = cell(pastframes,1);
    prev_bump_roi = [];
    bump_roi = [];
    y_shift = 0;
    x_shift = 0;
    usfac = 100; %tolerance for dft registration
    ev_ho_dict = zeros(21700,1); %at most 21600 frame long video, just add some extra frames, will be chopped (?)
    ev_ho_dict(curr_frame) = bgvidindex;
    
    %%% MAIN LOOP
    
    while hasFrame(v2) || ~isempty(movienames) %MAIN LOOP
        if ~hasFrame(v2) %change the video here
            disp('CHANGE VIDEO!');
            videoframe = 1;
            bgvidindex = bgvidindex+1;
            [ movienames, videoname ] = getnextvideo( movienames );
            v = VideoReader(videoname);%just for info
            [orig_background, clean_background, filledbg, level, abovethresh] = getbackground3(v,outer_boundary_mask,outer_boundary_line,level,0,abovethresh);
            bg_struct(bgvidindex).videoname = videoname;
            bg_struct(bgvidindex).videoframe = videoframe;
            bg_struct(bgvidindex).curr_frame = curr_frame;
            bg_struct(bgvidindex).region_rounded = region_rounded;
            bg_struct(bgvidindex).outer_boundary_mask = outer_boundary_mask;
            bg_struct(bgvidindex).outer_boundary_line = outer_boundary_line;
            bg_struct(bgvidindex).orig_background = orig_background;
            bg_struct(bgvidindex).clean_background = clean_background;
            bg_struct(bgvidindex).level = level;
            bg_struct(bgvidindex).abovethresh = abovethresh;
            [ ev_ho_x, ev_ho_y, lawn_limit_line] = geteventhorizon_rough(th,orig_background,lawn_limit_line); %event horizon is in whole frame coordinates
            ev_ho_crp_rel = [ev_ho_x-region_rounded(1) ev_ho_y-region_rounded(2)]; %event horizon in cropped coordinates
            lawn_limit_mask = poly2mask(lawn_limit_line(:,1),lawn_limit_line(:,2),size(filledbg,1),size(filledbg,2));
            fullyblurred_background = imgaussfilt(regionfill(filledbg,lawn_limit_mask),4);
            bg_struct(bgvidindex).fullyblurred_background = fullyblurred_background;
            bg_struct(bgvidindex).lawn_limit_mask_wf = lawn_limit_mask;
            bg_struct(bgvidindex).lawn_limit_line = lawn_limit_line;
            bg_struct(bgvidindex).ev_ho = [ev_ho_x ev_ho_y];
            bg_struct(bgvidindex).ev_ho_crp_rel = ev_ho_crp_rel;  % in cropped frame coordinates
            clear v2;
            v2 = VideoReader(videoname);%for tracking
            v2.CurrentTime = 0;
        end
        
        disp(['curr frame: ' num2str(curr_frame)]);
        frame_ui8 = rgb2gray(readFrame(v2)); %read in the next frame, uint8
        
        %%%% CHECK FOR VIDEO TRANSLATIONAL SHIFTS I.E. FROM PLATE BUMPS
        prev_bump_roi = bump_roi;
        bump_roi = imcrop(frame_ui8,crop_pos); %use this for checking if plate moved
        change_x = 0; change_y = 0;
        if curr_frame>1
            [output, ~] = dftregistration(fft2(prev_bump_roi),fft2(bump_roi),usfac); %register successive images
            add_y = output(3); add_x = output(4);
            if abs(add_y)>=1 %only take it seriously if the bump is at least 1 pixel
                y_shift = round(y_shift+add_y);
                change_y = 1;
            else
                change_y = 0;
            end
            if abs(add_x)>=1
                x_shift = round(x_shift+add_x);
                change_x = 1;
            else
                change_x = 0;
            end
        end
        %SHIFT EVERYTHING IF NECESSARY: OUTER BOUNDARY LINES, MASKS, BACKGROUND, EVENT
        %HORIZON (all masks and lines are defined with respect to the frame
        %size)
        if change_x || change_y
            region_rounded = [region_rounded(1)-x_shift region_rounded(2)-y_shift region_rounded(3)-abs(x_shift) region_rounded(4)-abs(y_shift)]; %the outer boundary line will always shrink if the plate is bumped --therefore minus the abs val
            tmpax = axes(); h = imellipse(tmpax,region_rounded);
            outer_boundary_line = getVertices(h); %this can move or change shape
            outer_boundary_mask = poly2mask(outer_boundary_line(:,1),outer_boundary_line(:,2),size(frame,1),size(frame,2));
            outer_boundary_line_crp_rel = [outer_boundary_line(:,1)-region_rounded(1) outer_boundary_line(:,2)-region_rounded(2)];%shift to the correct coordinates
            lawn_limit_line = [lawn_limit_line(:,1)-x_shift lawn_limit_line(:,2)-y_shift];
            lawn_limit_mask = poly2mask(lawn_limit_line(:,1),lawn_limit_line(:,2),size(frame,1),size(frame,2));
            ev_ho_x = ev_ho_x-x_shift; ev_ho_y = ev_ho_y-y_shift;
            ev_ho_crp_rel = [ev_ho_x-region_rounded(1) ev_ho_y-region_rounded(2)]; %event horizon in cropped coordinates, this should basically do nothing since you are adding and subtracting. just added for sake of clarity and to be sure everything matches up
            close all;
            %change these, save changes
            bgvidindex = bgvidindex+1;
            bg_struct(bgvidindex).videoname = videoname;
            bg_struct(bgvidindex).videoframe = videoframe;
            bg_struct(bgvidindex).curr_frame = curr_frame;
            bg_struct(bgvidindex).region_rounded = region_rounded;
            bg_struct(bgvidindex).outer_boundary_mask = outer_boundary_mask;
            bg_struct(bgvidindex).outer_boundary_line = outer_boundary_line;
            bg_struct(bgvidindex).orig_background = orig_background;
            bg_struct(bgvidindex).clean_background = clean_background;
            bg_struct(bgvidindex).level = level;
            bg_struct(bgvidindex).fullyblurred_background = fullyblurred_background;
            bg_struct(bgvidindex).abovethresh = abovethresh;
            bg_struct(bgvidindex).lawn_limit_mask_wf = lawn_limit_mask;
            bg_struct(bgvidindex).lawn_limit_line = lawn_limit_line;
            bg_struct(bgvidindex).ev_ho = [ev_ho_x ev_ho_y];
            bg_struct(bgvidindex).ev_ho_crp_rel = ev_ho_crp_rel;
        end
        %%%%
        
        frame = imcomplement(im2double(frame_ui8));
        %mask out abovethresh pixels from the background
        frame_ui8 = regionfill(frame_ui8,abovethresh.*outer_boundary_mask);
        frame = regionfill(frame,abovethresh.*outer_boundary_mask);
        %%%%%%%%%
        masked = frame.*outer_boundary_mask;
        cropped = imcrop(masked,region_rounded);
        bgsub = (imcomplement(cropped - imcrop(clean_background.*outer_boundary_mask,region_rounded)));
        thresh1 = imcomplement(im2bw(bgsub,level));
        thresh_cleaned1 = bwareaopen(thresh1,250);
        BWfinal1 = thresh_cleanup( thresh_cleaned1, 29, 3 );
        %%% Detect worms: %%%
        [area, centroids, bboxes] = step(hblob,BWfinal1); %image coordinates
        bboxes = double(bboxes);
        %%%%%%%%%%%%%%%%%%%%%
        %get a cleaned, worm-removed version of this image for extracting the
        %grayscale value (proxy for bacterial thickness)
        cropped_gs = get_clean_grayscale( frame_ui8, region_rounded, outer_boundary_mask, BWfinal1, bboxes, fullyblurred_background );
        
        %%% ASSIGN GRAYSCALE VALUE ASSOCIATED WITH THE HEAD AND TAIL
        lookback = 0; %boolean -- whether to look back or not
        %this is a first in first out queue
        if curr_frame<=pastframes
            gs_lookback{curr_frame} = cropped_gs;
            pt1_g = NaN;
            pt2_g = NaN;
            lookback = 0;
        else %start looking back
            %figure out which frame to use to get the grayscale head value
            cropped_gs_detect = cell(size(centroids,1),1);
            for k = 1:size(centroids,1) %get one past frame for each centroid detected
                cent_lookback = getprevframe4gs(tracks,tracksThatLeft,curr_frame,centroids(k,:),pastframes,cropped);
                numlookback = pastframes-cent_lookback+1;
                cropped_gs_detect{k} = gs_lookback{numlookback}; %1 is pastframes ago
            end
            gs_lookback(1) = [];%remove the first element
            gs_lookback{pastframes} = cropped_gs;%add the new element
            lookback = 1;
        end
        %%%%%%%%%
        
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
            [~,dist,~] = distance2curve(outer_boundary_line_crp_rel,centroids(j,:));
            if lookback
                curr_lookback_frame = cropped_gs_detect{j};
            end
            pt1_g = NaN; %default values, only to be overridden if lookback is active
            pt2_g = NaN;
            if dist > 40 %allowed to get spline at this distance!
                %using Ev's code for segmentation
                try
                    [worm,vWorm,errNum,errMsg] = segWorm_Elias(wormcrop,wormcrop_orig,curr_frame, 0, 1);
                    if isempty(errNum) %everything's groovy
                        pt1 = worm.skeleton.pixels(1,:);
                        pt1 = [pt1(2)+bbox_crop(1) pt1(1)+bbox_crop(2)]; %flip x and y to go from image to plotting
                        pt2 = worm.skeleton.pixels(end,:);
                        pt2 = [pt2(2)+bbox_crop(1) pt2(1)+bbox_crop(2)];
                        if lookback
                            pt1_g = curr_lookback_frame(pt1(2),pt1(1)); %get the grayscale value associated with the head and tail - remember this is image coordinates, have to reverse again!
                            pt2_g = curr_lookback_frame(pt2(2),pt2(1));
                        end
                    elseif ~isempty(worm) && errNum~= 101 && errNum ~= 102 && errNum ~= 103 && errNum ~= 105 %an error was raised but it's still fine
                        pt1 = worm.skeleton.pixels(1,:);
                        pt1 = [pt1(2)+bbox_crop(1) pt1(1)+bbox_crop(2)];
                        pt2 = worm.skeleton.pixels(end,:);
                        pt2 = [pt2(2)+bbox_crop(1) pt2(1)+bbox_crop(2)];
                        if lookback
                            pt1_g = curr_lookback_frame(pt1(2),pt1(1)); %get the grayscale value associated with the head and tail - remember this is image coordinates, have to reverse again!
                            pt2_g = curr_lookback_frame(pt2(2),pt2(1));
                        end
                    elseif isempty(worm) && ~isempty(vWorm.skeleton.pixels) %only save the endpoints, these can be helpful
                        worm = Inf;
                        pt1 = vWorm.skeleton.pixels(1,:);
                        pt1 = [pt1(2)+bbox_crop(1) pt1(1)+bbox_crop(2)];
                        pt2 = vWorm.skeleton.pixels(end,:);
                        pt2 = [pt2(2)+bbox_crop(1) pt2(1)+bbox_crop(2)];
                        if lookback
                            pt1_g = curr_lookback_frame(pt1(2),pt1(1)); %get the grayscale value associated with the head and tail - remember this is image coordinates, have to reverse again!
                            pt2_g = curr_lookback_frame(pt2(2),pt2(1));
                        end
                    else %the worm could not be segmented at all
                        worm = Inf;
                        pt1 = [NaN NaN];
                        pt2 = [NaN NaN];
                        %pt1_g, pt2_g are already specified with the default
                        %values above
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
            ev_ho_dict(curr_frame) = bgvidindex;
            continue;
        end
        %%% MAIN TRACKING FUNCTIONS
        tracks = predictNewLocationsOfTracks_wholevideo(tracks);
        
        [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment_wholevideo(tracks,centroids);
        
        if length(unassignedTracks) >= 1
            warning('WE LOST THE TRACK (NOT OUT OF BOUNDS)');
            disp(curr_frame);
        end
        
        tracks = updateAssignedTracks_wholevideo(tracks, assignments, x_shift, y_shift, centroids, bboxes, cropworms, cropworms_orig, splines, end1s, end2s, end1s_g, end2s_g, worms, curvatures, posture_angles, curr_frame, videoframe, videoname);
        
        [tracks, tracksThatLeft] = updateUnassignedTracks_wholevideo(tracks,tracksThatLeft,unassignedTracks);
        
        [tracks, nextId] = createNewTracks_wholevideo(tracks,nextId,x_shift,y_shift,unassignedDetections,centroids,bboxes,cropworms,cropworms_orig,splines,end1s, end2s, end1s_g, end2s_g, worms,curvatures, posture_angles, curr_frame,videoframe,videoname);

        curr_frame = curr_frame + 1;
        videoframe = videoframe +1;
        ev_ho_dict(curr_frame) = bgvidindex;
    end %END MAIN LOOP
    if curr_frame<length(ev_ho_dict) %chop it down to the correct size
        ev_ho_dict = ev_ho_dict(1:curr_frame);
    end
    
    %%% POST PROCESSING
    allTracks = [tracksThatLeft,tracks];
    WorthSavingTracksThatLeft_idx = [allTracks(:).age]'>9; %double check that all tracks are still long enough -- sometimes the last track is too short to be processed
    allTracks = allTracks(WorthSavingTracksThatLeft_idx);
    [TRACKS, EXIT_STRUCT, POKE_STRUCT] = tracks_postprocessing2( allTracks, ev_ho_dict, bg_struct, pixpermm );
    
    %%% SAVING!
    fields2remove = {'cropworm','cropworm_orig','worm'}; %<--- USE THIS ONE FOR ALL FAST DATA MANIPULATIONS
    TRACKS_slim = rmfield(TRACKS,fields2remove);
    
    save([videoname '_' lawn_string '_BACKGROUND.mat'],'bg_struct');
    save([videoname '_' lawn_string '_FINAL' '.mat'],'TRACKS_slim','EXIT_STRUCT','POKE_STRUCT','bg_struct');
    save([videoname '_' lawn_string '_TRACKS_CROPWORMS.mat'],'TRACKS','-v7.3');
    
end
end
