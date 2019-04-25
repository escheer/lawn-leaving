function track_QC_071218(stat_int, before_sec, after_sec)
%This function opens each .mat file in the folder and decides whether to
%include it in subsequent analyses based on the following criteria:
% 1. there was a worm in the video
% 2. the worm was visible for at least half of the video
% 3. the worm gets into the lawn within the first 20 minutes of the video
% 4. there aren't too many tracks after removing garbage tracks
% 5. the plate wasn't bumped more than 10 pixels in either direction during
% the video

%what's new in this version:
%assume you have bg_struct saved in _FINAL.mat!
%now call get_head_pokes2.m which calculates the headpokes+reversals


% for at least a cumulative minute's worth of time
% if these criteria are met, we copy the .mat to a subfolder called "passQC"

timenow = datestr(now,'mm_dd_yy');

[~, pathname, ~] = uigetfile({'*'});
cd(pathname);
files = dir('*FINAL.mat');

% out=regexp(pathname,'\','split');
% foldername = out{end-1};

foldername = '';
mkdir([foldername '_passQC']);
subdir = [pwd '\' foldername '_passQC'];

for i = 1:length(files)
    
    keep = true; %boolean whether to keep this file
    resave = false;
    
    disp(files(i).name);
    %load everythiing -- this takes a long time
%     tracks = load(files(i).name);
%     TRACKS_slim = tracks.TRACKS_slim;
%     EXIT_STRUCT = tracks.EXIT_STRUCT;
%     POKE_STRUCT = tracks.POKE_STRUCT;
%     bg_struct = tracks.bg_struct;

    %just load TRACKS_slim
    TRACKS_slim = load(files(i).name,'TRACKS_slim');
    TRACKS_slim = TRACKS_slim.TRACKS_slim; %always a struct within a struct


    
    %check if there are any tracks that need to be removed -- such as a
    %piece of dust, etc. this check should also be performed by
    %tracks_postprocessing2.m but I will put it here for now.
    bbox = {TRACKS_slim.bbox}';
    area = cell2mat(cellfun(@(x) mean(x(:,3).*x(:,4),1),bbox,'UniformOutput',false)); %get the average area of the bounding box for the worm for each track
    
    area_cutoff = 750;
    numtracks_cutoff = 150;
    if sum(area>area_cutoff)~=length(area)
        TRACKS_slim = TRACKS_slim(area>area_cutoff);
        
        if length(TRACKS_slim)>numtracks_cutoff %if you still have a ton of tracks after removing obvious crap, the whole thing is probably crap
            keep = false;
            disp('TOO MANY TRACKS! DISCARD!');
        else
            bg_struct = load(files(i).name,'bg_struct');
            bg_struct = bg_struct.bg_struct;
            %%%%%%%%% recalculate lawn entry and exit, remake plots
            [ TRACKS_slim, FRAMES_IN_LAWN, ENTER_FRAMES, ENTER_TRKS, EXIT_FRAMES, EXIT_TRKS, EXIT_INTS, ALIGNED_EXIT_IND, EXIT_TRK_KEY, INTS_IN_BY_TRACK, INTS_OUT_BY_TRACK,INTS_IN_GLOBAL,INTS_OUT_GLOBAL, IN_INT_TRK_KEY, OUT_INT_TRK_KEY] = get_enter_exit_events2( TRACKS_slim, before_sec, after_sec ); %get lawn entries and exits and intervals for each track surrounding the lawn exit event
            EXIT_STRUCT.FRAMES_IN_LAWN = FRAMES_IN_LAWN;
            EXIT_STRUCT.ENTER_FRAMES = ENTER_FRAMES;
            EXIT_STRUCT.ENTER_TRKS = ENTER_TRKS;
            EXIT_STRUCT.EXIT_FRAMES = EXIT_FRAMES;
            EXIT_STRUCT.EXIT_TRKS = EXIT_TRKS;
            EXIT_STRUCT.EXIT_INTS = EXIT_INTS;
            EXIT_STRUCT.ALIGNED_EXIT_IND = ALIGNED_EXIT_IND;
            EXIT_STRUCT.EXIT_TRK_KEY = EXIT_TRK_KEY;
            EXIT_STRUCT.INTS_IN_BY_TRACK = INTS_IN_BY_TRACK;
            EXIT_STRUCT.INTS_OUT_BY_TRACK = INTS_OUT_BY_TRACK;
            EXIT_STRUCT.INTS_IN_GLOBAL = INTS_IN_GLOBAL;
            EXIT_STRUCT.INTS_OUT_GLOBAL = INTS_OUT_GLOBAL;
            EXIT_STRUCT.IN_INT_TRK_KEY = IN_INT_TRK_KEY;
            EXIT_STRUCT.OUT_INT_TRK_KEY = OUT_INT_TRK_KEY;
            
            [TRACKS_slim, POKE_INTS_GLOBAL, POKE_INTS_BY_TRACK, POKE_TRK_KEY, POKE_PEAK_IDX_GLOBAL, POKE_PEAK_IDX_BY_TRACK, POKE_DIST_MINSUBTRACT, POKE_IS_REV, POKE_RAD_DIST, POKE_EVHO_DIST, POKE_GS] = get_head_pokes2( TRACKS_slim, bg_struct, INTS_OUT_BY_TRACK, OUT_INT_TRK_KEY );
            POKE_STRUCT.POKE_INTS_GLOBAL = POKE_INTS_GLOBAL;
            POKE_STRUCT.POKE_INTS_BY_TRACK = POKE_INTS_BY_TRACK;
            POKE_STRUCT.POKE_TRK_KEY = POKE_TRK_KEY;
            POKE_STRUCT.POKE_PEAK_IDX_GLOBAL = POKE_PEAK_IDX_GLOBAL;
            POKE_STRUCT.POKE_DIST_MINSUBTRACT = POKE_DIST_MINSUBTRACT;
            POKE_STRUCT.POKE_IS_REV = POKE_IS_REV;
            POKE_STRUCT.POKE_PEAK_IDX_BY_TRACK = POKE_PEAK_IDX_BY_TRACK;
            POKE_STRUCT.POKE_RAD_DIST = POKE_RAD_DIST;
            POKE_STRUCT.POKE_EVHO_DIST = POKE_EVHO_DIST;
            POKE_STRUCT.POKE_GS = POKE_GS;
            
%             startmin = 20; endmin = 60;
            NUMWORMS = 1;
%             stat_int = [startmin*60*3 endmin*60*3]; %frames between 20 minutes and 60 minutes since the start of the assay
            pixpermm = 112;
            
            foldername = [files(i).name '_QCedit'];
            curr_dir = pwd;
            
            %plot LAWN LEAVING over time
            [EXITS_DURING_INTERVAL, EXIT_COUNT_OVERTIME, EXIT_RATE_STATIC, ENTERS_DURING_INTERVAL, ENTER_COUNT_OVERTIME, ENTER_RATE_STATIC] = plot_LL_overtime2( ENTER_FRAMES, EXIT_FRAMES, NUMWORMS, FRAMES_IN_LAWN, stat_int, foldername, true, curr_dir);
            EXIT_STRUCT.EXITS_DURING_INTERVAL = EXITS_DURING_INTERVAL;
            EXIT_STRUCT.EXIT_COUNT_OVERTIME = EXIT_COUNT_OVERTIME;
            EXIT_STRUCT.EXIT_RATE_STATIC = EXIT_RATE_STATIC;
            EXIT_STRUCT.ENTERS_DURING_INTERVAL = ENTERS_DURING_INTERVAL;
            EXIT_STRUCT.ENTER_COUNT_OVERTIME = ENTER_COUNT_OVERTIME;
            EXIT_STRUCT.ENTER_RATE_STATIC = ENTER_RATE_STATIC;
            
            %plot HEAD POKES over time (now with a new category : HP+REV)
            [HEADPOKES_ALL, POKE_COUNT_OVERTIME_ALL, POKE_RATE_STATIC_ALL, AVG_POKE_DIST_ALL, HPREV, POKE_COUNT_OVERTIME_HPREV, POKE_RATE_STATIC_HPREV, AVG_POKE_DIST_HPREV] = plot_HP_overtime2( POKE_PEAK_IDX_GLOBAL, POKE_IS_REV, POKE_DIST_MINSUBTRACT, NUMWORMS, FRAMES_IN_LAWN, stat_int, pixpermm, foldername, true, curr_dir );
            POKE_STRUCT.HEADPOKES_DURING_INTERVAL_ALL = HEADPOKES_ALL;
            POKE_STRUCT.POKE_COUNT_OVERTIME_ALL = POKE_COUNT_OVERTIME_ALL;
            POKE_STRUCT.POKE_RATE_STATIC_ALL = POKE_RATE_STATIC_ALL;
            POKE_STRUCT.AVG_POKE_DIST_ALL = AVG_POKE_DIST_ALL;
            POKE_STRUCT.HEADPOKES_DURING_INTERVAL_REV = HPREV;
            POKE_STRUCT.POKE_COUNT_OVERTIME_REV = POKE_COUNT_OVERTIME_HPREV;
            POKE_STRUCT.POKE_RATE_STATIC_REV = POKE_RATE_STATIC_HPREV;
            POKE_STRUCT.AVG_POKE_DIST_ALL_REV = AVG_POKE_DIST_HPREV;
            %%%%%%%%%
            resave = true;
        end
    end
    
    if isempty(TRACKS_slim) %if there are no tracks, don't keep it
        keep = false;
    else %otherwise check other ways
        totalVisible = 0;
        totalInLawn = 0;
        totalInLawnbefore20 = 0;
        totalFrames = 10800; %disregard frames after 1 hour for the two hour videos I used to take
        for j = 1:length(TRACKS_slim)
            fa = TRACKS_slim(j).framesActive;
            fa_1 = fa(1);
            fa_end = fa(end);
            if fa_1 > totalFrames
                continue;
            end
            
            totalVisible = totalVisible + TRACKS_slim(j).totalVisibleCount;
            
            startidx = 1;
            stopidx = length(TRACKS_slim(j).centroidinlawn);
            if fa_1 <= totalFrames && fa_end > totalFrames
                stopidx = find(fa==10800);
                totalVisible = totalVisible - (fa_end-totalFrames); %remove overflow
            end
            totalInLawn = totalInLawn + nansum(TRACKS_slim(j).centroidinlawn(startidx:stopidx));
            
            if fa_1<3600 %20 minutes
                if fa_end > 3600
                    stopidx20 = find(fa==3600);
                else
                    stopidx20 = length(TRACKS_slim(j).centroidinlawn);
                end
                totalInLawnbefore20 = totalInLawnbefore20 + nansum(TRACKS_slim(j).centroidinlawn(startidx:stopidx20));
            end
        end
        %MUST BE...
        if totalVisible/totalFrames<0.50 % 50% visible within the hour
            disp('was not visible at least 50%.');
            keep = false;
        end
        
        %         if totalInLawn/totalFrames<0.50 % 50% in lawn during the hour
        %             disp('was not in lawn at least 50%.');
        %             keep = false;
        %         end
        
        if totalInLawnbefore20<180      % in the lawn for at least a minute within first 20 minutes
            disp('was not in lawn for at least 1 minute before 20 minutes.');
            keep = false;
        end
        
        %throw out any data that came from a video that was bumped more
        %than 10 pixels in either direction -- this data can be unreliable.
        mvmnt_cutoff = 10;
        xshift = {TRACKS_slim.x_shift}';
        xshiftmorethanx = cellfun(@(x) sum(x>mvmnt_cutoff)>0,xshift,'UniformOutput',true);
        yshift = {TRACKS_slim.y_shift}';
        yshiftmorethanx = cellfun(@(x) sum(x>mvmnt_cutoff)>0,yshift,'UniformOutput',true);
        if sum(xshiftmorethanx)>0 || sum(yshiftmorethanx)>0
            disp('plate was bumped more than 10 pixels in either direction!');
            keep = false;
        end
        
    end
    
    if keep
        if ~resave
%             copyfile(files(i).name,subdir);

%             cd(subdir);
%             save(files(i).name,'TRACKS_slim','EXIT_STRUCT','POKE_STRUCT','bg_struct');
%             cd(pathname);

            command = ['copy ' files(i).name ' ' subdir];
            %for long filenames, this may work better
%             command = ['robocopy /E ' pwd ' ' subdir ' ' files(i).name];
            status = dos(command);
            if status ~= 0
                disp('not copied successfully!');
            end
            
        else
            bg_struct = load(files(i).name,'bg_struct');
            bg_struct = bg_struct.bg_struct;
            cd(subdir);
            save([files(i).name(1:end-25) '_QC_EDIT_' timenow '_FINAL.mat'],'TRACKS_slim','EXIT_STRUCT','POKE_STRUCT','bg_struct'); %just save the same tracks that you loaded -- which includes head_pokes and leaving_events
            cd(pathname);
        end
        
    else
        disp('discard!');
    end
    
%     pause(0.1); %this may help not clog up page-ins?
end

end

