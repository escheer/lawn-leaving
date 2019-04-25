function REANALYZE_BATCH_010219(subfoldersToReanalyze)
warnid = 'images:initSize:adjustingMag';
warning('off',warnid);

timenow = datestr(now,'mm_dd_yy');
lawn_string = ['010219_REGEN_BATCH_' timenow];
p=10; %use the 10th percentile to generate the background
num_gaussian = 10; %default was 3

camera_params_path = which('tracking_params\camera_params_02_02_2019.mat');%locate the tracking parameters file
load(camera_params_path, 'camera_params'); %load in data

%if pixpermm is not supplied, get this value from camera_params

files = dir('*.avi');
uncutmovies = files(not([files.isdir]'));
moviename = uncutmovies(1).name;
cameracode = moviename(1:2);
camera_names = {'s1','s2','s3','s4','s5','s6','c1','c2','c3','c4','c5','c6'};
recognized = strcmp(camera_names,cameracode);
while ~sum(recognized)
    disp('Camera number not recognized.');
    cameracode = input('PLEASE SUPPLY CAMERA CODE (e.g. ''c1'')');
    recognized = strcmp(camera_names,cameracode);
end
switch cameracode
    case 's1'
        level = camera_params.s1.level;
        pixpermm = camera_params.s1.pixpermm;
    case 's2'
        level = camera_params.s2.level;
        pixpermm = camera_params.s2.pixpermm;
    case 's3'
        level = camera_params.s3.level;
        pixpermm = camera_params.s3.pixpermm;
    case 's4'
        level = camera_params.s4.level;
        pixpermm = camera_params.s4.pixpermm;
    case 's5'
        level = camera_params.s5.level;
        pixpermm = camera_params.s5.pixpermm;
    case 's6'
        level = camera_params.s6.level;
        pixpermm = camera_params.s6.pixpermm;
    case 'c1'
        level = camera_params.c1.level;
        pixpermm = camera_params.c1.pixpermm;
    case 'c2'
        level = camera_params.c2.level;
        pixpermm = camera_params.c2.pixpermm;
    case 'c3'
        level = camera_params.c3.level;
        pixpermm = camera_params.c3.pixpermm;
    case 'c4'
        level = camera_params.c4.level;
        pixpermm = camera_params.c4.pixpermm;
    case 'c5'
        level = camera_params.c5.level;
        pixpermm = camera_params.c5.pixpermm;
    case 'c6'
        level = camera_params.c6.level;
        pixpermm = camera_params.c6.pixpermm;
    otherwise
        disp('Camera data cannot be found!');
        level = input('Input Threshold Level (0-1)');
        pixpermm = input('Input Number of Pixels per MM');
end


[filename, outer_directory, ~] = uigetfile({'*CONFIG.mat'});
cd(outer_directory);
batch_struct = struct();
load(filename,'batch_struct');


%% TRACKING!
for vid = 1:length(batch_struct)
    tmp = batch_struct(vid).path;
    split = strsplit(tmp,'\');
    pathname = [outer_directory split{end} '\'];
    %go to the correct directory
    cd(pathname);
    disp('------------------------------------------------------');
    disp('------------------------------------------------------');
    disp(['NOW TRACKING: ' pathname]);
    disp('------------------------------------------------------');
    disp('------------------------------------------------------');
    %% LOAD IN BACKGROUND AND TRACKS_CROPWORMS.mat FILE
    localfiles = dir('*');
    localfilenames = {localfiles.name}';
    
    folders = localfiles([localfiles.isdir]'); %get the subdirectories
    folders = folders(~contains({folders.name}','.'));
    
    if exist('subfoldersToReanalyze','var') %if you give it a list of videos to track, just track those.
        if ~sum(strcmp(split{end},subfoldersToReanalyze))>0
            disp('Skip!');
            continue;
        end
    end
    
    %check if this file has already been reanalyzed, if so then skip it!
    if ~isempty(folders)&&length(folders)==1
        cd(folders.name);
        tmp = dir('*_REGEN_*.mat');
        if ~isempty(tmp)
            disp('Already Re-analyzed! Skip!');
            continue;
        end
    end
    
    BackgroundIdx = contains(localfilenames,'_BACKGROUND.mat');
    BackgroundFinished = sum(BackgroundIdx)>0;
    TracksCropwormsIdx = contains(localfilenames,'_TRACKS_CROPWORMS.mat');
    TracksCropwormsFinished = sum(TracksCropwormsIdx)>0;
    
    if ~BackgroundFinished || ~TracksCropwormsFinished
        error('In order to re-analyze old data, you must have already generated a BACKGROUND.mat and TRACKS_CROPWORMS.mat file!');
    end
    %Initialize variables that will be loaded.
    bg_struct = struct();
    TRACKS = struct();
    
    bgFile = localfilenames{BackgroundIdx};
    disp('Loading bg_struct...');
    load(bgFile,'bg_struct');
    tcFile = localfilenames{TracksCropwormsIdx};
    disp('Loading TRACKS_CROPWORMS (this will take a minute)...');
    load(tcFile,'TRACKS');
    
    %% STEP 1, RE-CALCULATE THE EVENT HORIZON FROM WITHIN THE LAWN-LIMIT LINE
    % AND OTHER TRACKING PARAMETERS
    origbgstruct = bg_struct; %copy original bg_struct. bg_struct will be edited.
    %Initialize variables used in shared functions.
    videoname = [];
    croppedcleanedBg = [];
    
    %% STEP 2, RE-ANALYZE THESE TRACKS BASED ON NEW TRACKING INFORMATION.
    NEW_TRACKS = TRACKS; %copy TRACKS, these fields will be edited.
    curr_bgvidindex = 0; %impossible value, will have to be updated on first loop.
    for trk = 1:length(NEW_TRACKS)
        disp(['TRACKNUM = ' num2str(trk)]);
        disp('---------------------------');
        track = NEW_TRACKS(trk);
        for idx = 1:track.age
            disp(['TRACKIDX = ' num2str(idx)]);
            if track.bgvidindex(idx) ~= curr_bgvidindex %switch tracking parameters when bgvidindex changes!
                curr_bgvidindex = track.bgvidindex(idx);
                updateBgStruct(curr_bgvidindex);
            end
            % Read in original cropworms (not thresholded), subtract newly
            % calculated backgrounds, threshold, proceed with tracking steps...
            bbox = track.bbox(idx,:);
            bbox_crop = [bbox(1)-5 bbox(2)-5 bbox(3)+10 bbox(4)+10];
            cw_orig = imcomplement(track.cropworm_orig{idx}); %you have to invert it!
            local_bg = imcrop(croppedcleanedBg, bbox_crop);
            bgsub = (imcomplement(cw_orig - local_bg));
            thresh = imcomplement(im2bw(bgsub,level));
            
            %since we are thresholding this worm freshly, let's pad the image
            %with 5 pixels of zeros on each side (this will avoid any clearborder
            %errors
            thresh = padarray(thresh,[5 5]);
            x_offset = bbox_crop(1)-5;
            y_offset = bbox_crop(2)-5;
            
            pixremoval_scale_factor = (500/96);
            thresh_cleaned = bwareaopen(thresh,round(pixremoval_scale_factor*pixpermm));
            
            wormcrop = thresh_cleanup_110818(thresh_cleaned); %this just erodes and dilates, no contour smoothing.
            if sum(sum(wormcrop))==0 %if the worm image ends up empty, don't try to segment it
                continue;
            end
            try
                [wormskel, pt1, pt2, omega] = segment_worm();
            catch
                continue; %if for some reason the segmentation fails on this frame, just skip it.
            end
            
            if ~isnan(wormskel) %if the worm has been successfully segmented, get it's spline and associated values
                [spline, curvature, posture_angle] = getWormSpline010219(wormskel);
            else
                spline = NaN(49,2);
                curvature = NaN(49,1);
                posture_angle = NaN(1,48); %the only oddball that got transposed. can fix later if desired.
            end
            
            %edit these fields in NEW_TRACKS.
            track.cropworm{idx} = wormcrop;
            track.end1(idx,:) = pt1;
            track.end2(idx,:) = pt2;
            track.omega(idx) = omega;
            track.spline{idx} = spline;
            track.curvature{idx} = curvature;
            track.posture_angle{idx} = posture_angle;
        end
        NEW_TRACKS(trk) = track; %update this track.
    end
    
    %% STEP 3. RE-PROCESS THESE TRACKS. SAVE NEW INFORMATION.
    %     mkdir([videoname(1:end-4) '_' lawn_string]); %save in a
    %     subdirectory. (can make filenames too long)
    %     cd([videoname(1:end-4) '_' lawn_string]);
    mkdir(lawn_string);%save in a subdirectory.
    cd(lawn_string);
    
    [TRACKS, EXIT_STRUCT, POKE_STRUCT, SUMMARY_STRUCT] = tracks_postprocessing5( NEW_TRACKS, [], bg_struct, pixpermm, [videoname(1:end-4) '_' lawn_string] );
    
    fields2remove = {'cropworm','cropworm_orig'};
    TRACKS_slim = rmfield(TRACKS,fields2remove);
    
    save([videoname '_' lawn_string '_BACKGROUND.mat'],'bg_struct');
    save([videoname '_' lawn_string '_FINAL' '.mat'],'TRACKS_slim','EXIT_STRUCT','POKE_STRUCT','SUMMARY_STRUCT','bg_struct');
    save([videoname '_' lawn_string '_SUMMARY.mat'],'SUMMARY_STRUCT');
    save([videoname '_' lawn_string '_TRACKS_CROPWORMS.mat'],'TRACKS','-v7.3');
    cd(pathname);
    
end

%NESTED FUNCTIONS
    function updateBgStruct(bgvidindex)
        videoname = origbgstruct(bgvidindex).videoname;
        v = VideoReader([pathname videoname]);
        number_of_frames = v.NumberOfFrames;
        last_frame = imadjust(rgb2gray(im2double(read(v,number_of_frames))));
        
        %get outer boundary mask, outer boundary line, region rounded, outer_boundary_line_crp_rel from old bg_struct
        outer_boundary_mask = origbgstruct(bgvidindex).outer_boundary_mask;
        outer_boundary_line = origbgstruct(bgvidindex).outer_boundary_line;
        region_rounded = origbgstruct(bgvidindex).region_rounded;
        
        %re-determine the background
        [orig_background, clean_background] = getbackground7_prctile_cleanup(v, level, p, num_gaussian, pixpermm, outer_boundary_mask);
        
        %re-generate the event horizon, taking into account the lawn-limit line as
        %an additional constraint (this avoids the need for user input).
        lawn_limit_line = origbgstruct(bgvidindex).lawn_limit_line;
        %         [ ev_ho_x, ev_ho_y, lawn_limit_line] = geteventhorizon010819(orig_background,outer_boundary_mask, lawn_limit_line);
        [ ev_ho_x, ev_ho_y, lawn_limit_line] = geteventhorizon022119(orig_background, outer_boundary_mask, lawn_limit_line, pixpermm);
        ev_ho_crp_rel = [ev_ho_x-region_rounded(1) ev_ho_y-region_rounded(2)]; %event horizon in cropped coordinates
        lawn_limit_mask_wf = poly2mask(lawn_limit_line(:,1),lawn_limit_line(:,2),size(last_frame,1),size(last_frame,2));
        fullyblurred_background = imgaussfilt(regionfill(orig_background,lawn_limit_mask_wf),4);
        
        %to be used locally for tracking.
        croppedcleanedBg = imcrop(clean_background.*outer_boundary_mask,region_rounded);
        
        %update information in config bg_struct
        bg_struct(bgvidindex).videoname = videoname;
        bg_struct(bgvidindex).videoframe = origbgstruct(bgvidindex).videoframe;
        bg_struct(bgvidindex).curr_frame = origbgstruct(bgvidindex).curr_frame;
        bg_struct(bgvidindex).region_rounded = region_rounded;
        bg_struct(bgvidindex).outer_boundary_mask = outer_boundary_mask;
        bg_struct(bgvidindex).outer_boundary_line = outer_boundary_line;
        bg_struct(bgvidindex).orig_background = orig_background;
        bg_struct(bgvidindex).clean_background = clean_background;
        bg_struct(bgvidindex).level = level;
        bg_struct(bgvidindex).fullyblurred_background = fullyblurred_background;
        bg_struct(bgvidindex).lawn_limit_mask_wf = lawn_limit_mask_wf;  % in whole frame coordinates
        bg_struct(bgvidindex).lawn_limit_line = lawn_limit_line;     % in whole frame coordinates
        bg_struct(bgvidindex).ev_ho = [ev_ho_x ev_ho_y];             % in whole frame coordinates
        bg_struct(bgvidindex).ev_ho_crp_rel = ev_ho_crp_rel;         % in cropped frame coordinates
    end

    function [wormskel, pt1, pt2, omega] = segment_worm()
        [wormskel, pt1, pt2, omega] = segWormOnly_Elias2(wormcrop, x_offset, y_offset, pixpermm);
    end


end
