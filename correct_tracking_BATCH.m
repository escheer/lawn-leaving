function correct_tracking_BATCH(NUMWORMS, startmin, endmin)
clc; close all;
%CORRECT_TRACKING_BATCH.m This function corrects mistakes in tracking and
%analyzed lawn leaving events and head pokes for all subfolders of the
%selected directory.


[~, pathname, ~] = uigetfile({'*'}); %this is nice because you can see which ones are in progress
cd(pathname);
files = dir('*');
folders = files([files.isdir]'); %get the subdirectories
folders = folders(~contains({folders.name}','.')); %get rid of the extra folders starting with '.' and '..'
outer_directory = pwd;
str = outer_directory ;
idx = strfind(str,'\') ;
foldername = str(idx(end)+1:end) ;

for i = 1:length(folders) %go into each subdirectory and get all the necessary information for tracking
    cd(outer_directory);
    curr_dir = folders(i).name;
    cd(curr_dir);
    bglist = dir('*_BACKGROUND.mat');
    bgname = {bglist.name}';
    if length(bgname)==0 %if missing file or multiple files, skip folder
        continue;
    end
    [~,I] = max([bglist(:).datenum]); %get the newest file
    bg_struct = load(bgname{I});
    bg_struct = bg_struct.bg_struct; %BACKGROUND STRUCT
    
    trklist = dir('*_SLIM.mat');
    if length(trklist)==0 %if no files here, skip
        continue;
    end
    trkname = {trklist.name}';
    [~,I] = max([trklist(:).datenum]); %get the newest file
    TRACKS = load(trkname{I}); %TRACKS BEFORE CORRECTION AND ANALYSIS
    TRACKS = TRACKS.allTracks_slim;
    
    %LOOK FOR TRACKS OF NON-WORMS and REMOVE THEM. - USE MEAN SQUARED
    %DISPLACEMENT
    msds = {TRACKS.msd}';
    avg_msd = cellfun(@nanmean,msds);
%     remove = 0;
%     badtracks = avg_msd<=1e-5;
%     if sum(badtracks>0)
%         plotting_head_grayscale_on_background_ALLTRACKS(TRACKS,bg_struct(1).region_rounded(1),bg_struct(1).region_rounded(2),bg_struct(1).orig_background,.75,1.1);
%         remove = input('ARE THERE DUST TRACKS? (1/0)');
%         close all;
%     end
    
%     if remove
        TRACKS = TRACKS(avg_msd>1e-5); %remove any tracks with inordinately low MSD
%     end
    
    %REPLACE THE SMOOTHED CENTROID, HEAD, TAIL WITH A LESS SMOOTHED VERSION (WINDOW = 3),
    %FIX DERIVED VARIABLES
    for i = 1:length(TRACKS)
        centroid = TRACKS(i).centroid;
        centroid_smooth = [movmean(centroid(:,1),3,'omitnan') movmean(centroid(:,2),3,'omitnan')];
        TRACKS(i).centroid_smooth = centroid_smooth;
        
        pixpermm = 112;
        [msd, msd_smth] = get_msd( centroid_smooth, pixpermm, 3 );
        TRACKS(i).msd = msd;
        TRACKS(i).msd_smooth = msd_smth;
        
        %GET SPEED AND ANGULAR SPEED
        set1 = centroid_smooth(1:end-3,:);
        set2 = centroid_smooth(4:end,:); %1 second ahead for calculation
        speed = [zeros(3,1); diag(pdist2(set1,set2))./pixpermm]; % D is the displacement between successive centroid positions in millimeters
        TRACKS(i).speed = speed; % (mm/sec) absolute value of the velocity (scale this per 1/3 second = 3fps)
        [angspeed,~] = getAngularSpeed2(centroid_smooth,3);
        TRACKS(i).angspeed = angspeed;
        
        % FIX ROAMING DWELLING ANNOTATIONS
        trans = [0.967, 0.033; 0.043, 0.957];
        emis = [0.89, 0.11; 0.30, 0.70];
        cutoff = 1250;
        binSize = 5*3;%5 seconds = 15 frames
        [~, expStates, estTR, estE, actualratio] = getHMMStates(TRACKS(i),binSize,TRACKS(i).centroidinlawn,cutoff,trans,emis);
        TRACKS(i).roamdwell = expStates;
        TRACKS(i).actualratio = actualratio;
        TRACKS(i).estTR = estTR;
        TRACKS(i).estE = estE;
    end
    
    TRACKS = fix_head_tail( TRACKS, bg_struct ); %correct mistakes in head and tail annotation, this also re-smoothes head, tail, speed, rearranges gs values to match head/tail fixes, calls get forward and reverse again, calls countBlobsInOut again
    [ TRACKS, FRAMES_IN_LAWN, ENTER_FRAMES, ENTER_TRKS, EXIT_FRAMES, EXIT_TRKS, EXIT_INTS, ALIGNED_EXIT_IND, EXIT_TRK_KEY, INTS_IN_BY_TRACK, INTS_OUT_BY_TRACK,INTS_IN_GLOBAL,INTS_OUT_GLOBAL, IN_INT_TRK_KEY, OUT_INT_TRK_KEY] = get_enter_exit_events2( TRACKS, 45, 15 ); %get lawn entries and exits and intervals for each track surrounding the lawn exit event
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
    
    [TRACKS, POKE_INTS_GLOBAL, POKE_INTS_BY_TRACK, POKE_TRK_KEY, POKE_PEAK_IDX_GLOBAL, POKE_PEAK_IDX_BY_TRACK, POKE_DIST_MINSUBTRACT, POKE_RAD_DIST, POKE_EVHO_DIST, POKE_GS] = get_head_pokes( TRACKS, bg_struct, 0.970, INTS_OUT_BY_TRACK, OUT_INT_TRK_KEY ); %get head pokes
    POKE_STRUCT.POKE_INTS_GLOBAL = POKE_INTS_GLOBAL;
    POKE_STRUCT.POKE_INTS_BY_TRACK = POKE_INTS_BY_TRACK;
    POKE_STRUCT.POKE_TRK_KEY = POKE_TRK_KEY;
    POKE_STRUCT.POKE_PEAK_IDX_GLOBAL = POKE_PEAK_IDX_GLOBAL;
    POKE_STRUCT.POKE_DIST_MINSUBTRACT = POKE_DIST_MINSUBTRACT;
    POKE_STRUCT.POKE_PEAK_IDX_BY_TRACK = POKE_PEAK_IDX_BY_TRACK;
    POKE_STRUCT.POKE_RAD_DIST = POKE_RAD_DIST;
    POKE_STRUCT.POKE_EVHO_DIST = POKE_EVHO_DIST;
    POKE_STRUCT.POKE_GS = POKE_GS;
    
    stat_int = [startmin*60*3 endmin*60*3]; %frames between 20 minutes and 60 minutes since the start of the assay

    %plot LAWN LEAVING over time
    [EXITS_DURING_INTERVAL, EXIT_COUNT_OVERTIME, EXIT_RATE_STATIC, ENTERS_DURING_INTERVAL, ENTER_COUNT_OVERTIME, ENTER_RATE_STATIC] = plot_LL_overtime( ENTER_FRAMES, EXIT_FRAMES, NUMWORMS, FRAMES_IN_LAWN, stat_int, foldername);
    EXIT_STRUCT.EXITS_DURING_INTERVAL = EXITS_DURING_INTERVAL;
    EXIT_STRUCT.EXIT_COUNT_OVERTIME = EXIT_COUNT_OVERTIME;
    EXIT_STRUCT.EXIT_RATE_STATIC = EXIT_RATE_STATIC;
    EXIT_STRUCT.ENTERS_DURING_INTERVAL = ENTERS_DURING_INTERVAL;
    EXIT_STRUCT.ENTER_COUNT_OVERTIME = ENTER_COUNT_OVERTIME;
    EXIT_STRUCT.ENTER_RATE_STATIC = ENTER_RATE_STATIC;
    
    %plot HEAD POKES over time
    [HEADPOKES_DURING_INTERVAL, POKE_COUNT_OVERTIME, POKE_RATE_STATIC, AVG_POKE_DIST] = plot_HP_overtime( POKE_PEAK_IDX_GLOBAL, POKE_DIST_MINSUBTRACT, NUMWORMS, FRAMES_IN_LAWN, stat_int, pixpermm, foldername );
    POKE_STRUCT.HEADPOKES_DURING_INTERVAL = HEADPOKES_DURING_INTERVAL;
    POKE_STRUCT.POKE_COUNT_OVERTIME = POKE_COUNT_OVERTIME;
    POKE_STRUCT.POKE_RATE_STATIC = POKE_RATE_STATIC;
    POKE_STRUCT.AVG_POKE_DIST = AVG_POKE_DIST;
    
    close all;
    cd(outer_directory); %go back to outer directory to save the final .mats -- makes them easier to find
    TRACKS_FILENAME = trkname{I};
    save([TRACKS_FILENAME(1:end-15) '_FINAL.mat'],'TRACKS', 'EXIT_STRUCT', 'POKE_STRUCT');
end


end