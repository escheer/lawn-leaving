function correct_tracking_BATCH_newHPs_043018(NUMWORMS, startmin, endmin, out_dir, plotflag)
clc; close all;
%CORRECT_TRACKING_BATCH.m This function corrects mistakes in tracking and
%analyzed lawn leaving events and head pokes for all subfolders of the
%selected directory.
pixpermm = 112;

[~, pathname, ~] = uigetfile({'*'}); %this is nice because you can see which ones are in progress
cd(pathname);

% files = dir('*');
% folders = files([files.isdir]'); %get the subdirectories
% folders = folders(~contains({folders.name}','.')); %get rid of the extra folders starting with '.' and '..'
% outer_directory = pwd;
% str = outer_directory ;
% idx = strfind(str,'\') ;
% trackname = str(idx(end)+1:end) ;

% for i = 1:length(folders) %go into each subdirectory and get all the necessary information for tracking
%     cd(outer_directory);
%     curr_dir = folders(i).name;
%     cd(curr_dir);
%     bglist = dir('*_BACKGROUND.mat');
%     bgname = {bglist.name}';
%     if length(bgname)==0 %if missing file or multiple files, skip folder
%         continue;
%     end
%     [~,I] = max([bglist(:).datenum]); %get the newest file
%     bg_struct = load(bgname{I});
%     bg_struct = bg_struct.bg_struct; %BACKGROUND STRUCT
%
%     trklist = dir('*_FINAL.mat');
%     if length(trklist)==0 %if no files here, skip
%         continue;
%     end

trklist = dir('*_FINAL.mat');
trkname = {trklist.name}';

for i = 1:length(trklist)
    cd(pathname);
    TRACKS_everything = load(trkname{i}); %TRACKS BEFORE CORRECTION AND ANALYSIS
    disp(trkname{i});
    %     if i == 10
    %         disp('debug');
    %     end
    if isfield(TRACKS_everything,'TRACKS')
        TRACKS = TRACKS_everything.TRACKS;
    elseif isfield(TRACKS_everything,'TRACKS_slim')
        TRACKS = TRACKS_everything.TRACKS_slim;
    else
        error('must have either TRACKS or TRACKS_slim!');
    end
    
    try
        bg_struct = TRACKS_everything.bg_struct;
    catch
    end
    %     (since we have already run this once before)
    %check that roaming and dwelling data are there:
    if ~isfield(TRACKS,'roamdwell_hmm')
        msds = {TRACKS.msd}';
        avg_msd = cellfun(@nanmean,msds);
        TRACKS = TRACKS(avg_msd>1e-5); %remove any tracks with inordinately low MSD
        for trkidx = 1:length(TRACKS)
            %SMOOTH CENTROID, GET MEAN SQUARED DISPLACEMENT
            centroid = TRACKS(trkidx).centroid;
            centroid_smooth = [movmean(centroid(:,1),3,'omitnan') movmean(centroid(:,2),3,'omitnan')];
            TRACKS(trkidx).centroid_smooth = centroid_smooth;
            
            [msd, msd_smth] = get_msd( centroid_smooth, pixpermm, 3 );
            TRACKS(trkidx).msd = msd;
            TRACKS(trkidx).msd_smooth = msd_smth;
            
            %GET SPEED AND ANGULAR SPEED
            set1 = centroid_smooth(1:end-3,:);
            set2 = centroid_smooth(4:end,:); %1 second ahead for calculation
            speed = [zeros(3,1); diag(pdist2(set1,set2))./pixpermm]; % D is the displacement between successive centroid positions in millimeters
            TRACKS(trkidx).speed = speed; % (mm/sec) absolute value of the velocity (scale this per 1/3 second = 3fps)
            [angspeed,~] = getAngularSpeed3(centroid_smooth,1);
            TRACKS(trkidx).angspeed = angspeed;
            
            %my numbers from N2 OD2
            trans = [0.9645    0.0355; 0.0802    0.9198];
            emis =  [0.9790    0.0210; 0.5448    0.4552];
            cutoff = 35;
            x_offset = 2.5;
            binSize = 10*3;%10 seconds = 30 frames
            train = false;
            
            if TRACKS(trkidx).age>=binSize %if the track length is longer than binSize, categories roaming and dwelling states
                [expSeq, expStates, ~, ~, actualratio] = getHMMStates3(TRACKS(trkidx),binSize,TRACKS(trkidx).centroidinlawn,cutoff,trans,emis,x_offset,train);
                TRACKS(trkidx).roamdwell_2d = expSeq;
                TRACKS(trkidx).roamdwell_hmm = expStates;
                TRACKS(trkidx).actualratio = actualratio;
            else
                TRACKS(trkidx).roamdwell_2d = NaN(TRACKS(trkidx).age,1);
                TRACKS(trkidx).roamdwell_hmm = NaN(TRACKS(trkidx).age,1);
                TRACKS(trkidx).actualratio = NaN(TRACKS(trkidx).age,1);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    before_sec = 300; %5 minutes before
    after_sec = 60; %1 minute after
    [ TRACKS, FRAMES_IN_LAWN, ENTER_FRAMES, ENTER_TRKS, EXIT_FRAMES, EXIT_TRKS, EXIT_INTS, ALIGNED_EXIT_IND, EXIT_TRK_KEY, INTS_IN_BY_TRACK, INTS_OUT_BY_TRACK,INTS_IN_GLOBAL,INTS_OUT_GLOBAL, IN_INT_TRK_KEY, OUT_INT_TRK_KEY] = get_enter_exit_events2( TRACKS, before_sec, after_sec ); %get lawn entries and exits and intervals for each track surrounding the lawn exit event
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
    
    [ TRACKS, POKE_INTS_GLOBAL, POKE_INTS_BY_TRACK, POKE_TRK_KEY, POKE_PEAK_IDX_GLOBAL, POKE_PEAK_IDX_BY_TRACK, POKE_DIST_MINSUBTRACT, POKE_IS_REV, POKE_RAD_DIST, POKE_EVHO_DIST, POKE_GS] = get_head_pokes2_forcorrection( TRACKS, INTS_OUT_BY_TRACK, OUT_INT_TRK_KEY );
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
    
    %     startmin = 20; endmin = 60;
%     pixpermm = 112;
    %     NUMWORMS = 1;
    stat_int = [startmin*60*3 endmin*60*3]; %frames between 20 minutes and 60 minutes since the start of the assay
    %     curr_dir = pwd;
    %     cd('..');
    %     str = pwd ;
    %     idx = strfind(str,'\') ;
    %     foldername = str(idx(end)+1:end) ;
    %     cd(curr_dir);
    
    str = trkname{i};
    expression = 'od[125]';
    idx_OD = regexp(str,expression);
    expression = '_c[123456]_';
    idx_cnum = regexp(str,expression);
    if idx_OD>idx_cnum
        idx = idx_OD;
    elseif idx_cnum> idx_OD
        idx = idx_cnum;
    else
        error('either Ruby or Elias naming convention must be operative!');
    end
    trackname = str(1:idx(1)+27);
    
    %plot LAWN LEAVING over time
    [EXITS_DURING_INTERVAL, EXIT_COUNT_OVERTIME, EXIT_RATE_STATIC, ENTERS_DURING_INTERVAL, ENTER_COUNT_OVERTIME, ENTER_RATE_STATIC] = plot_LL_overtime2( ENTER_FRAMES, EXIT_FRAMES, NUMWORMS, FRAMES_IN_LAWN, stat_int, trackname, false, out_dir);
    EXIT_STRUCT.EXITS_DURING_INTERVAL = EXITS_DURING_INTERVAL;
    EXIT_STRUCT.EXIT_COUNT_OVERTIME = EXIT_COUNT_OVERTIME;
    EXIT_STRUCT.EXIT_RATE_STATIC = EXIT_RATE_STATIC;
    EXIT_STRUCT.ENTERS_DURING_INTERVAL = ENTERS_DURING_INTERVAL;
    EXIT_STRUCT.ENTER_COUNT_OVERTIME = ENTER_COUNT_OVERTIME;
    EXIT_STRUCT.ENTER_RATE_STATIC = ENTER_RATE_STATIC;
    
    %plot HEAD POKES over time (now with a new category : HP+REV)
    [HEADPOKES_ALL, POKE_COUNT_OVERTIME_ALL, POKE_RATE_STATIC_ALL, AVG_POKE_DIST_ALL, HPREV, POKE_COUNT_OVERTIME_HPREV, POKE_RATE_STATIC_HPREV, AVG_POKE_DIST_HPREV] = plot_HP_overtime2( POKE_PEAK_IDX_GLOBAL, POKE_IS_REV, POKE_DIST_MINSUBTRACT, NUMWORMS, FRAMES_IN_LAWN, stat_int, pixpermm, trackname, false, out_dir );
    POKE_STRUCT.HEADPOKES_DURING_INTERVAL_ALL = HEADPOKES_ALL;
    POKE_STRUCT.POKE_COUNT_OVERTIME_ALL = POKE_COUNT_OVERTIME_ALL;
    POKE_STRUCT.POKE_RATE_STATIC_ALL = POKE_RATE_STATIC_ALL;
    POKE_STRUCT.AVG_POKE_DIST_ALL = AVG_POKE_DIST_ALL;
    POKE_STRUCT.HEADPOKES_DURING_INTERVAL_REV = HPREV;
    POKE_STRUCT.POKE_COUNT_OVERTIME_REV = POKE_COUNT_OVERTIME_HPREV;
    POKE_STRUCT.POKE_RATE_STATIC_REV = POKE_RATE_STATIC_HPREV;
    POKE_STRUCT.AVG_POKE_DIST_ALL_REV = AVG_POKE_DIST_HPREV;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    close all;
    cd(out_dir); %go to the out_dir for saving
    if exist('bg_struct')
        save([trackname '_mod04302018_FINAL.mat'],'TRACKS', 'EXIT_STRUCT', 'POKE_STRUCT','bg_struct');
    else
        save([trackname '_mod04302018_nobgstruct_FINAL.mat'],'TRACKS', 'EXIT_STRUCT', 'POKE_STRUCT');
    end
end


end