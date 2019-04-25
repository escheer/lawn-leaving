function updateV1tracks_022018()
%This function re-runs getAngularSpeed3, getHMMStates3 and
%get_enter_exit_events2 to keep old analysis consistent with new analysis
%including roaming and dwelling and keeping short lawn leaving events.

timenow = datestr(now,'mm_dd_yy');
lawn_string = ['022018_UPDATED_' timenow];

[~, pathname, ~] = uigetfile({'*'});
cd(pathname);
files = dir('*FINAL.mat');

mkdir('020218_updated');
subdir = [pwd '\' '020218_updated'];

for i = 1:length(files)
    disp(files(i).name);
    tracks = load(files(i).name);
    TRACKS = tracks.TRACKS;
    EXIT_STRUCT = tracks.EXIT_STRUCT;
    POKE_STRUCT = tracks.POKE_STRUCT;
    
    %remove some old ones
    if isfield(TRACKS, 'roamdwell')
        TRACKS = rmfield(TRACKS,'roamdwell');
    end
    if isfield(TRACKS, 'estTR')
        TRACKS = rmfield(TRACKS,'estTR');
    end
    if isfield(TRACKS, 'estE')
        TRACKS = rmfield(TRACKS,'estE');
    end
    
    for j = 1:length(TRACKS)
        disp(j)
        %fix off by ones on track length
        age = size(TRACKS(j).centroid,1);
        if TRACKS(j).age < age
            age = TRACKS(j).age;
        end
        %also check for tracks that end by skipping to a new location far from
        %where they just were -- these final observations should be cut off.
        A = TRACKS(j).centroid(1:age,:);
        dist = sqrt( sum( abs( diff( A ) ).^2, 2 ) );
        if dist(end)>mean(dist(1:end-1))+2*std(dist(1:end-1)) %if the track moves too far at the last frame, cut it back
            age = age-1;
        end
        
        TRACKS(j).id = TRACKS(j).id;
        TRACKS(j).bgvidindex = TRACKS(j).bgvidindex(1:age);
        TRACKS(j).x_shift = TRACKS(j).x_shift(1:age);
        TRACKS(j).y_shift = TRACKS(j).y_shift(1:age);
        
        TRACKS(j).bbox = TRACKS(j).bbox(1:age,:);
        
%         TRACKS(j).end1 = TRACKS(j).end1(1:age,:); %nice to keep track of these in case you need them later
%         TRACKS(j).end2 = TRACKS(j).end2(1:age,:);
%         TRACKS(j).end1_g = TRACKS(j).end1_g(1:age);
%         TRACKS(j).end2_g = TRACKS(j).end2_g(1:age);
        
        TRACKS(j).spline = TRACKS(j).spline(1:age);
        TRACKS(j).curvature = TRACKS(j).curvature(1:age);
        TRACKS(j).posture_angle = TRACKS(j).posture_angle(1:age);
        
        TRACKS(j).head = TRACKS(j).head(1:age,:);
        TRACKS(j).head_smooth = TRACKS(j).head_smooth(1:age,:);
        TRACKS(j).tail = TRACKS(j).tail(1:age,:);
        TRACKS(j).tail_smooth = TRACKS(j).tail_smooth(1:age,:);
        TRACKS(j).centroid = TRACKS(j).centroid(1:age,:);
        TRACKS(j).centroid_smooth = TRACKS(j).centroid_smooth(1:age,:);
        
        TRACKS(j).speed = TRACKS(j).speed(1:age);
        TRACKS(j).speed_smooth = TRACKS(j).speed_smooth(1:age);
        TRACKS(j).msd = TRACKS(j).msd(1:age);
        TRACKS(j).msd_smooth = TRACKS(j).msd_smooth(1:age);
        
        TRACKS(j).head_gs = TRACKS(j).head_gs(1:age);
        TRACKS(j).tail_gs = TRACKS(j).tail_gs(1:age);
        
        TRACKS(j).age = age;
        TRACKS(j).consecutiveInvisibleCount = TRACKS(j).consecutiveInvisibleCount;
        TRACKS(j).totalVisibleCount = age;
        TRACKS(j).framesActive = TRACKS(j).framesActive(1:age);
        TRACKS(j).videoname = TRACKS(j).videoname(1:age);
        TRACKS(j).videoframe = TRACKS(j).videoframe(1:age);
        
        TRACKS(j).omega = TRACKS(j).omega(1:age);
        TRACKS(j).forward = TRACKS(j).forward(1:age);
        TRACKS(j).reverse = TRACKS(j).reverse(1:age);
        
        TRACKS(j).fullyinlawn = TRACKS(j).fullyinlawn(1:age);
        TRACKS(j).headinlawn = TRACKS(j).headinlawn(1:age);
        TRACKS(j).centroidinlawn = TRACKS(j).centroidinlawn(1:age);
        TRACKS(j).tailinlawn = TRACKS(j).tailinlawn(1:age);
        
        TRACKS(j).lawn_entries = TRACKS(j).lawn_entries(1:age);
        TRACKS(j).lawn_exits = TRACKS(j).lawn_exits(1:age);
        TRACKS(j).head_pokes = TRACKS(j).head_pokes(1:age);
        TRACKS(j).radial_dist = TRACKS(j).radial_dist(1:age);
        TRACKS(j).ev_ho_dist = TRACKS(j).ev_ho_dist(1:age);
        
        %fix angular speed
        centroid_smooth = TRACKS(j).centroid_smooth;
        [angspeed,~] = getAngularSpeed3(centroid_smooth,1);
        TRACKS(j).angspeed = angspeed; %update with new method
        
        %fix roaming and dwelling with N2 OD2 values from Feb 2018
        trans = [0.9645    0.0355; 0.0802    0.9198];
        emis =  [0.9790    0.0210; 0.5448    0.4552];
        cutoff = 35;
        x_offset = 2.5;
        binSize = 10*3;%10 seconds = 30 frames
        train = false;
        
        if TRACKS(j).age>=binSize %if the track length is longer than binSize, categories roaming and dwelling states
            [expSeq, expStates, ~, ~, actualratio] = getHMMStates3(TRACKS(j),binSize,TRACKS(j).centroidinlawn,cutoff,trans,emis,x_offset,train);
            TRACKS(j).roamdwell_2d = expSeq;
            TRACKS(j).roamdwell_hmm = expStates;
            TRACKS(j).actualratio = actualratio;
        else
            TRACKS(j).roamdwell_2d = NaN(TRACKS(j).age,1);
            TRACKS(j).roamdwell_hmm = NaN(TRACKS(j).age,1);
            TRACKS(j).actualratio = NaN(TRACKS(j).age,1);
        end
    end
    
    %fix enter and exit events
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
    
    %don't need to fix the head pokes
    
    startmin = 20; endmin = 60;
    NUMWORMS = 1;
    stat_int = [startmin*60*3 endmin*60*3]; %frames between 20 minutes and 60 minutes since the start of the assay

    foldername = files(i).name(1:end-25);
    
    cd(subdir); %change to subdirectory to save plots and new .mat files
    %plot LAWN LEAVING over time
    [EXITS_DURING_INTERVAL, EXIT_COUNT_OVERTIME, EXIT_RATE_STATIC, ENTERS_DURING_INTERVAL, ENTER_COUNT_OVERTIME, ENTER_RATE_STATIC] = plot_LL_overtime( EXIT_STRUCT.ENTER_FRAMES, EXIT_STRUCT.EXIT_FRAMES, NUMWORMS, EXIT_STRUCT.FRAMES_IN_LAWN, stat_int, foldername);
    EXIT_STRUCT.EXITS_DURING_INTERVAL = EXITS_DURING_INTERVAL;
    EXIT_STRUCT.EXIT_COUNT_OVERTIME = EXIT_COUNT_OVERTIME;
    EXIT_STRUCT.EXIT_RATE_STATIC = EXIT_RATE_STATIC;
    EXIT_STRUCT.ENTERS_DURING_INTERVAL = ENTERS_DURING_INTERVAL;
    EXIT_STRUCT.ENTER_COUNT_OVERTIME = ENTER_COUNT_OVERTIME;
    EXIT_STRUCT.ENTER_RATE_STATIC = ENTER_RATE_STATIC;
    
    %plot HEAD POKES over time.
    pixpermm = 112;
    [HEADPOKES_DURING_INTERVAL, POKE_COUNT_OVERTIME, POKE_RATE_STATIC, AVG_POKE_DIST] = plot_HP_overtime( POKE_STRUCT.POKE_PEAK_IDX_GLOBAL, POKE_STRUCT.POKE_DIST_MINSUBTRACT, NUMWORMS, EXIT_STRUCT.FRAMES_IN_LAWN, stat_int, pixpermm, foldername );
    POKE_STRUCT.HEADPOKES_DURING_INTERVAL = HEADPOKES_DURING_INTERVAL;
    POKE_STRUCT.POKE_COUNT_OVERTIME = POKE_COUNT_OVERTIME;
    POKE_STRUCT.POKE_RATE_STATIC = POKE_RATE_STATIC;
    POKE_STRUCT.AVG_POKE_DIST = AVG_POKE_DIST;
    
    %rename TRACKS TRACKS_slim to be consistent with other data
    TRACKS_slim = TRACKS;
    %save new .mat file
    save([foldername '_' lawn_string '_FINAL' '.mat'],'TRACKS_slim','EXIT_STRUCT','POKE_STRUCT'); %this file will be same as new ones except it does not include the bg_struct (oh well...)
    close all;
    cd(pathname); %return to outer folder to go on to the next file
end

end

