%THIS CODE RE-SEGMENTS PREVIOUSLY TRACKED VIDEOS
%Elias Scheer
%11-08-18. This version attempts to use the 10th percentile as the
%background instead of the mean.

function [TRACKS, EXIT_STRUCT, POKE_STRUCT, TRACKS_slim] = tracking_enter_exit_RESEGMENT_BATCH(recalc_bg)
warnid = 'images:initSize:adjustingMag';
warning('off',warnid);

timenow = datestr(now,'mm_dd_yy');
lawn_string = ['110818_RESEGMENTED_BATCH_' timenow];

level = 0.905;

p=10; %use the 10th percentile to generate the background

th = 0.001; %for event horizon detection

[filename, pathname, ~] = uigetfile({'*_CONFIG.mat'}); %this is nice because you can see which ones are in progress
cd(pathname);

T = load(filename,'batch_struct');%this loads batch_struct
batch_struct = T.('batch_struct');

for vid = 1:length(batch_struct) %outer loop, over each subfolder where there are cropped videos and tracking files
    tic;
    %go to the correct directory
    cd(batch_struct(vid).path);
    %%% GET DATA FROM BATCH STRUCT
    pixpermm = batch_struct(vid).pixpermm;
    %%% LOAD IN BG_STRUCT (this will duplicate some information from batch_struct -- but it has data for each sub-video)
    bg_file = ls('*_BACKGROUND.mat');
    if isempty(bg_file)
        error('This video has not been tracked! Cannot retrack!');
    elseif size(bg_file,1)>1 %if there is more than one background file, choose which one you want to use
        [bg_file, ~, ~] = uigetfile({'*_BACKGROUND.mat'});
    end
    disp('LOADING BACKGROUND FILE...');
    T = load(bg_file,'bg_struct');
    bg_struct_OLD = T.('bg_struct'); %rename bg_struct to bg_struct_OLD to free up the variable name
    
    %%% LOAD IN CROPWORM DATA FROM LAST TIME MOVIE WAS TRACKED
    cropworm_tracks_file = ls('*_TRACKS_CROPWORMS.mat');
    if isempty(cropworm_tracks_file)
        error('This video has not been tracked! Cannot retrack!');
    elseif size(cropworm_tracks_file,1)>1 %if there is more than one TRACKS file, choose which one you want to use
        [cropworm_tracks_file, ~, ~] = uigetfile({'*_TRACKS_CROPWORMS.mat'});
    end
    avi_char_idx = strfind(cropworm_tracks_file,'.avi')+3;
    disp('LOADING CROPWORMS FILE...');
    T = load(cropworm_tracks_file,'TRACKS'); %this will take awhile, rename it to have a copy of the original
    TRACKS_old = T.('TRACKS');
    
    %START RE-SEGMENTING
    bg_struct = bg_struct_OLD;% make a copy of bg_struct_OLD, fields will be replaced if necessary
    TRACKS = TRACKS_old; %make a copy TRACKS_old, these values will be updated after re-segmentation
    
    currvid = ''; %initialize this so a VideoReader is made on the first loop
    for i = 1:length(TRACKS) %loop over TRACKS
        disp(['TRACK ' num2str(i) ':']);
        track = TRACKS(i);
        videoname = track.videoname;
        videoframe = track.videoframe;
        %pre-allocate cell arrays for the fields that will be replaced directly
        cropworm = cell(track.age,1);
        cropworm_orig = cell(track.age,1);
        cropworm_normIllum = cell(track.age,1);
        reseg_bbox_crop = NaN(track.age,4); %coordinates of crop worm in whole frame coordinates
        spline = cell(track.age,1);
        curvature = cell(track.age,1);
        posture_angle = cell(track.age,1);
        worm = cell(track.age,1);
        %make a dictionary of pointsoutsideworm for each frame
        points_outside_worm = cell(track.age,1);
        %         cw = []; cw_orig = [];
        for j = 1:track.age %loop over indices of that track
%             if j==21
%                 disp('debug');
%             end
            curr_frame = videoframe(j);
            disp(['idx ' num2str(j) ' / videoframe ' num2str(curr_frame)]);
            vn = videoname{j};
            
            if ~strcmp(currvid,vn)  %change videos
                disp('CHANGING VIDEOS!');
                v = VideoReader(vn);
                if recalc_bg        %recalculate background and stuff related to that if specified
                    region_rounded = bg_struct_OLD(track.bgvidindex(j)).region_rounded;
                    lawn_limit_line = bg_struct_OLD(track.bgvidindex(j)).lawn_limit_line; %you need this so you don't have to manually redraw
                    disp('RE-CALCULATING BACKGROUND...');
                    recalculate_background();
                else
                    orig_background = bg_struct_OLD(track.bgvidindex(j)).orig_background;
                    clean_background = bg_struct_OLD(track.bgvidindex(j)).clean_background;
                    ev_ho = bg_struct_OLD(track.bgvidindex(j)).ev_ho;
                    ev_ho_crp_rel = bg_struct_OLD(track.bgvidindex(j)).ev_ho_crp_rel;
                    lawn_limit_line = bg_struct_OLD(track.bgvidindex(j)).lawn_limit_line;
                    lawn_limit_mask_wf = bg_struct_OLD(track.bgvidindex(j)).lawn_limit_mask_wf;
                    fullyblurred_background = bg_struct_OLD(track.bgvidindex(j)).fullyblurred_background;
                end
                %UPDATE BG_STRUCT (most fields remain the same)
                bg_struct(track.bgvidindex(j)).orig_background = orig_background;
                bg_struct(track.bgvidindex(j)).clean_background = clean_background;
                bg_struct(track.bgvidindex(j)).ev_ho = ev_ho;
                bg_struct(track.bgvidindex(j)).ev_ho_crp_rel = ev_ho_crp_rel;
                bg_struct(track.bgvidindex(j)).lawn_limit_line = lawn_limit_line;
                bg_struct(track.bgvidindex(j)).lawn_limit_mask_wf = lawn_limit_mask_wf;
                bg_struct(track.bgvidindex(j)).fullyblurred_background = fullyblurred_background;
                currvid = vn;%update
            end
            %Now, extract cropped worm images, subtract the background,
            %resegment, etc.
            frame_x_offset = bg_struct_OLD(track.bgvidindex(j)).region_rounded(1);
            frame_y_offset = bg_struct_OLD(track.bgvidindex(j)).region_rounded(2);
            bbox = track.bbox(j,:);
            bbox_crop = [bbox(1)-5 bbox(2)-5 bbox(3)+10 bbox(4)+10];
            bbox_crop_wf = [bbox_crop(1)+frame_x_offset bbox_crop(2)+frame_y_offset bbox_crop(3) bbox_crop(4)];%in wholeframe coordinates
            cw_orig = track.cropworm_orig{j};
            if ~(size(cw_orig,1)-bbox_crop(4)==1 && size(cw_orig,2)-bbox_crop(3)==1) %sometimes when the cropped image is at the edge of the well, this bug happens
                cropworm{j} = track.cropworm{j};
                cropworm_orig{j} = track.cropworm_orig{j};
                reseg_bbox_crop(j,:) = [NaN NaN NaN NaN];
                worm{j} = track.worm{j};
                spline{j} = track.spline{j};
                curvature{j} = track.curvature{j};
                posture_angle{j} = track.posture_angle{j};
                continue;
            end
           
            %perform background subtraction and thresholding
            bg_crop = imcrop(clean_background,bbox_crop_wf);
            bgsub = imcomplement(imcomplement(cw_orig) - bg_crop);
            thresh = imcomplement(im2bw(bgsub,level));
            thresh_cleaned = bwareaopen(thresh,250);
            cw = thresh_cleanup_RESEG(thresh_cleaned);
            cw = padarray(cw,[1 1],0,'both'); 
            %pad the cropped worm image to ensure the contour
            %never touches the boundary of the image. This increases the size of the image by
            %1 pixel in all directions (width+2, height+2)
            cw_orig = padarray(cw_orig,[1 1],'replicate','both'); %match it
            %this is written redundantly, but for clarity. 
            new_x_offset = bbox_crop_wf(1)-1; new_y_offset = bbox_crop_wf(2)-1; %compensate for these offsets
            bbox_crop_wf = [new_x_offset new_y_offset bbox_crop_wf(3)+2 bbox_crop_wf(4)+2]; %update this to match
            
            %get an illumination-normalized cropworm for extracting
            %grayscale
            invfilledbg = imcrop(imcomplement(fullyblurred_background),bbox_crop_wf);
            normIllumCropWorm = cw_orig./invfilledbg;
            
            %get the points in cropworm containing no worm (for calculating
            %head grayscale value)
            gsmask = imdilate(cw,strel('disk',7)); %expand further off the worm contour to define "points on worm"
            [pow_y, pow_x] = find(~gsmask);
            pointsoutsideworm = [pow_x pow_y]; %you can find out if a point is in this list by: find(pointsoutsideworm(:,1)==27 & pointsoutsideworm(:,2)==17)
            points_outside_worm{j} = [pointsoutsideworm(:,1)+new_x_offset pointsoutsideworm(:,2)+new_y_offset];

            %segment the worm and extract the grayscale value for the
            %endpoints
            [wrm, pt1, pt2, pt1_g, pt2_g] = segment_worm();

            if isstruct(wrm) %if the worm has been successfully segmented, get it's spline and associated values
                [spln, crvtr, pstr_angl] = getWormSpline2( wrm, new_x_offset, new_y_offset );
            else
                spln = NaN;
                crvtr = NaN;
                pstr_angl = NaN;
            end
            
            cropworm{j} = cw;
            cropworm_orig{j} = cw_orig;
            cropworm_normIllum{j} = normIllumCropWorm;
            reseg_bbox_crop(j,:) = bbox_crop_wf;
            worm{j} = wrm;
            spline{j} = spln;
            curvature{j} = crvtr;
            posture_angle{j} = pstr_angl;
            track.end1(j,:) = pt1;
            track.end2(j,:) = pt2;
            track.end1_g(j) = pt1_g;
            track.end2_g(j) = pt2_g;
        end
        TRACKS(i).cropworm = cropworm;
        TRACKS(i).cropworm_orig = cropworm_orig;
        TRACKS(i).cropworm_normIllum = cropworm_normIllum;
        TRACKS(i).reseg_bbox_crop = reseg_bbox_crop;
        TRACKS(i).worm = worm;
        TRACKS(i).spline = spline;
        TRACKS(i).curvature = curvature;
        TRACKS(i).posture_angle = posture_angle;
%         TRACKS(i) = track; %update this track
    end
    [TRACKS, EXIT_STRUCT, POKE_STRUCT] = tracks_postprocessing_RESEG( TRACKS, bg_struct, pixpermm );
    
    %%% SAVING!
    fields2remove = {'cropworm','cropworm_orig','worm'}; %<--- USE THIS ONE FOR ALL FAST DATA MANIPULATIONS
    TRACKS_slim = rmfield(TRACKS,fields2remove);
    
    save([cropworm_tracks_file(1:avi_char_idx) '_' lawn_string '_BACKGROUND.mat'],'bg_struct');
    save([cropworm_tracks_file(1:avi_char_idx) '_' lawn_string '_FINAL' '.mat'],'TRACKS_slim','EXIT_STRUCT','POKE_STRUCT','bg_struct');
    save([cropworm_tracks_file(1:avi_char_idx) '_' lawn_string '_TRACKS_CROPWORMS.mat'],'TRACKS','-v7.3');
    
    toc;
end


%%%% NESTED FUNCTIONS %%%%

    function recalculate_background()
        [orig_background, clean_background, filledbg, level] = getbackground4_prctile( v, level, p);
        [ ev_ho_x, ev_ho_y, lawn_limit_line] = geteventhorizon_rough(th,orig_background,lawn_limit_line);
        ev_ho = [ev_ho_x ev_ho_y];
        ev_ho_crp_rel = [ev_ho_x-region_rounded(1) ev_ho_y-region_rounded(2)]; %event horizon in cropped coordinates
        lawn_limit_mask_wf = poly2mask(lawn_limit_line(:,1),lawn_limit_line(:,2),size(orig_background,1),size(orig_background,2));
        fullyblurred_background = imgaussfilt(regionfill(filledbg,lawn_limit_mask_wf),4);
    end

    function [worm, pt1, pt2, pt1_g, pt2_g] = segment_worm()
        pt1_g = NaN; %default values, only to be overridden if a grayscale value can be found
        pt2_g = NaN;
        %using Ev's code
%         try
            [worm,vWorm,errNum,~] = segWorm_Elias(cw,cw_orig,curr_frame, 0, 1);
            if isempty(errNum) %everything's groovy
                pt1 = worm.skeleton.pixels(1,:);
                pt1 = [pt1(2) pt1(1)];%flip x and y
                pt1_g = get_grayscale(pt1);
                pt1 = [pt1(1)+new_x_offset pt1(2)+new_y_offset]; %add back offsets
                pt2 = worm.skeleton.pixels(end,:);
                pt2 = [pt2(2) pt2(1)];%flip x and y
                pt2_g = get_grayscale(pt2);
                pt2 = [pt2(1)+new_x_offset pt2(2)+new_y_offset]; %add back offsets
            elseif ~isempty(worm) && errNum~= 101 && errNum ~= 102 && errNum ~= 103 && errNum ~= 105 %an error was raised but it's still fine
                pt1 = worm.skeleton.pixels(1,:);
                pt1 = [pt1(2) pt1(1)];%flip x and y
                pt1_g = get_grayscale(pt1);
                pt1 = [pt1(1)+new_x_offset pt1(2)+new_y_offset]; %add back offsets
                pt2 = worm.skeleton.pixels(end,:);
                pt2 = [pt2(2) pt2(1)];%flip x and y
                pt2_g = get_grayscale(pt2);
                pt2 = [pt2(1)+new_x_offset pt2(2)+new_y_offset]; %add back offsets
            elseif isempty(worm) && ~isempty(vWorm.skeleton.pixels) %only save the endpoints, these can be helpful
                worm = Inf;
                pt1 = vWorm.skeleton.pixels(1,:);
                pt1 = [pt1(2) pt1(1)];%flip x and y
                pt1_g = get_grayscale(pt1);
                pt1 = [pt1(1)+new_x_offset pt1(2)+new_y_offset]; %add back offsets
                pt2 = vWorm.skeleton.pixels(end,:);
                pt2 = [pt2(2) pt2(1)];%flip x and y
                pt2_g = get_grayscale(pt2);
                pt2 = [pt2(1)+new_x_offset pt2(2)+new_y_offset]; %add back offsets
            else %the worm could not be segmented at all
                worm = Inf;
                pt1 = [NaN NaN];
                pt2 = [NaN NaN];
                %pt1_g, pt2_g are already specified with the default (NaN)
                %values above
            end
%         catch
%             warning('THERE WAS AN ERROR IN THE SEGMENTATION CODE!');
%             worm = Inf;
%             pt1 = [NaN NaN];
%             pt2 = [NaN NaN];
%         end
        
        function pt_g = get_grayscale(pt) 
            %get the grayscale value associated with each of the worm endpoints
            pt_g = NaN; %initialize this value; no way to look back on the first frame
            lookback_idx = j-1;%track index minus 1
            pointsfound = false;
            while ~pointsfound && lookback_idx>0 %look back through cropworm images to find a frame containing the endpoint without the worms body in it.
                prev_pow = points_outside_worm{lookback_idx};%wf coords
                %get the list of pixels containing the endpoint and its
                %neighbors
%                 [pt_and_neighbs,~] = get_pixel_neighbors(cw,pt); %pt must be in local indexing
%                 pt_and_neighbs = [pt_and_neighbs;pt]; %add back in the point itself
                pt_and_neighbs = pt;%just try with the point itself
                pt_and_neighbs_wf = [pt_and_neighbs(:,1)+new_x_offset pt_and_neighbs(:,2)+new_y_offset];%add back wf offsets to do the checking
                %check whether pixels associated with endpoints in the
                %current frame are points_outside_worm in previous frames.
                pf = zeros(size(pt_and_neighbs,1),1); %whether the point was found
                pf_val = NaN(size(pt_and_neighbs,1),1); %the grayscale value of the point if it was found
                isCurrHeadOutsidePrevWorm = sum(prev_pow(:,1)==pt(1)+new_x_offset & prev_pow(:,2)==pt(2)+new_y_offset);
                PrevImageBbox = reseg_bbox_crop(lookback_idx,:);
                % image vertices 1->2
                %                4<-3
                PrevImageVerts = [PrevImageBbox(1) PrevImageBbox(2); PrevImageBbox(1)+PrevImageBbox(3) PrevImageBbox(2);...
                                    PrevImageBbox(1)+PrevImageBbox(3) PrevImageBbox(2)+PrevImageBbox(4); PrevImageBbox(1) PrevImageBbox(2)+PrevImageBbox(4)];
                [in,on] = inpolygon(pt_and_neighbs_wf(:,1),pt_and_neighbs_wf(:,2),PrevImageVerts(:,1),PrevImageVerts(:,2));
                isCurrHeadInsidePrevImage = in | on;
                if isCurrHeadOutsidePrevWorm && isCurrHeadInsidePrevImage    %sum(isCurrHeadInsidePrevImage)>=3 %only look if at least 3 pixels of pt and neighbs is in the previous image
                    cw_normIllum_lkback = cropworm_normIllum{lookback_idx};
                    for n = 1:size(pt_and_neighbs,1)
                        PointIsInImage = pt_and_neighbs(n,2)<=size(cw_normIllum_lkback,1) && pt_and_neighbs(n,1)<=size(cw_normIllum_lkback,2);
                        PointIsOutsideWorm = sum(prev_pow(:,1)==pt_and_neighbs_wf(n,1) & prev_pow(:,2)==pt_and_neighbs_wf(n,2)); %will be 1 if the point was found, 0 if not
                        pf(n) = PointIsInImage && PointIsOutsideWorm;
                        
                        if (PointIsInImage && PointIsOutsideWorm) %if this point is outside the old worm and contained with the image, take its value
                            pf_val(n) = cw_normIllum_lkback(pt_and_neighbs(n,2),pt_and_neighbs(n,1));
                        end
                    end
                    if mean(pf)>=(1/3) %if at least 1/3 of the points+neighbors were found, take these points' grayscale values
                        pt_g = nanmean(pf_val);
                        pointsfound = true;
%                         disp(['grayscale = ' num2str(pt_g)]);
%                         h=overlay_line_on_image(cw_normIllum_lkback,pt_and_neighbs,10,'r'); pause(); close(h);
                    end
                end
                lookback_idx = lookback_idx-1;
            end
            
        end
    end
end