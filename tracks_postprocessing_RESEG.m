function [TRACKS, EXIT_STRUCT, POKE_STRUCT] = tracks_postprocessing_RESEG( allTracks, bg_struct, pixpermm )
%TRACKS_POSTPROCESSING.m This function takes in the current tracks and the
%tracks that finished throughout the video and extracts a few more
%behavioral parameters from them.

TRACKS = struct(...
    'id', {}, ...
    'bgvidindex',{},...
    'x_shift',{},...
    'y_shift',{},...
    'centroid', {}, ...
    'bbox', {}, ...
    'cropworm',{}, ...
    'cropworm_orig',{},...
    'cropworm_normIllum',{},...
    'reseg_bbox_crop',{},...
    'spline',{},...
    'curvature',{},...
    'posture_angle',{},...
    'worm',{},...
    'head',{},...
    'head_smooth',{},...
    'tail',{},...
    'tail_smooth',{},...
    'end1',{},...
    'end2',{},...
    'centroid_smooth',{},...
    'speed',{},...
    'speed_smooth',{},...
    'msd',{},...
    'msd_smooth',{},...
    'angspeed',{},...
    'actualratio',{},...
    'roamdwell_hmm',{},...
    'roamdwell_2d',{},...
    'head_gs',{},...
    'tail_gs',{},...
    'age', {}, ...
    'totalVisibleCount', {}, ...
    'consecutiveInvisibleCount', {}, ...
    'framesActive', {},...
    'videoname',{},...
    'videoframe',{},...
    'omega',{},...
    'forward',{},...
    'reverse',{},...
    'fullyinlawn',{},...
    'headinlawn',{},...
    'centroidinlawn',{},...
    'tailinlawn',{});

for i = 1:length(allTracks)
    disp(['curr track is ' num2str(i)]);
    % occasionally you end up with different number of centroid positions
    % than the age of the track and this should be fixed in the main
    % tracking code, but also check for it here and fix it.
    age = size(allTracks(i).centroid,1);
    if allTracks(i).age < age
        age = allTracks(i).age;
    end
    %also check for tracks that end by skipping to a new location far from
    %where they just were -- these final observations should be cut off.
    A = allTracks(i).centroid(1:age,:);
    dist = sqrt( sum( abs( diff( A ) ).^2, 2 ) );
    if dist(end)>mean(dist(1:end-1))+2*std(dist(1:end-1)) %if the track moves too far at the last frame, cut it back
        age = age-1;
    end
    %%% copy most fields directly to TRACKS
    TRACKS(i).id = allTracks(i).id;
    TRACKS(i).bgvidindex = allTracks(i).bgvidindex;
    TRACKS(i).x_shift = allTracks(i).x_shift(1:age);
    TRACKS(i).y_shift = allTracks(i).y_shift(1:age);
    TRACKS(i).centroid = allTracks(i).centroid(1:age,:);
    TRACKS(i).bbox = allTracks(i).bbox(1:age,:);
    
    TRACKS(i).cropworm = allTracks(i).cropworm(1:age);
    TRACKS(i).cropworm_orig = allTracks(i).cropworm_orig(1:age);
    TRACKS(i).cropworm_normIllum = allTracks(i).cropworm_normIllum(1:age);
    TRACKS(i).reseg_bbox_crop = allTracks(i).reseg_bbox_crop(1:age);%new field added -- these are the coordinates of the upper left hand corner of the new cropworm
    TRACKS(i).end1 = allTracks(i).end1(1:age,:); %nice to keep track of these in case you need them later
    TRACKS(i).end2 = allTracks(i).end2(1:age,:);
    TRACKS(i).end1_g = allTracks(i).end1_g(1:age);
    TRACKS(i).end2_g = allTracks(i).end2_g(1:age);
    TRACKS(i).worm = allTracks(i).worm(1:age);
    
    TRACKS(i).spline = allTracks(i).spline(1:age);
    TRACKS(i).curvature = allTracks(i).curvature(1:age);
    TRACKS(i).posture_angle = allTracks(i).posture_angle(1:age);
    TRACKS(i).age = age;
    TRACKS(i).consecutiveInvisibleCount = allTracks(i).consecutiveInvisibleCount;
    TRACKS(i).totalVisibleCount = age;
    TRACKS(i).framesActive = allTracks(i).framesActive(1:age);
    TRACKS(i).videoname = allTracks(i).videoname(1:age);
    TRACKS(i).videoframe = allTracks(i).videoframe(1:age);
    
    frames = TRACKS(i).framesActive;
%     bgvidindex = ev_ho_dict(frames); %GET RID OF THIS FOR RESEG FUNCTION.
%     TRACKS(i).bgvidindex = bgvidindex(1:age); %this can be used to index indto the bg_struct to access any of that information on a per frame basis
    %%%
    %SMOOTH CENTROID, GET MEAN SQUARED DISPLACEMENT
    centroid = TRACKS(i).centroid;
    centroid_smooth = [movmean(centroid(:,1),3,'omitnan') movmean(centroid(:,2),3,'omitnan')];
    TRACKS(i).centroid_smooth = centroid_smooth;
    
    [msd, msd_smth] = get_msd( centroid_smooth, pixpermm, 3 );
    TRACKS(i).msd = msd;
    TRACKS(i).msd_smooth = msd_smth;
    
    %GET SPEED AND ANGULAR SPEED
    set1 = centroid_smooth(1:end-3,:);
    set2 = centroid_smooth(4:end,:); %1 second ahead for calculation
    speed = [zeros(3,1); diag(pdist2(set1,set2))./pixpermm]; % D is the displacement between successive centroid positions in millimeters
    TRACKS(i).speed = speed; % (mm/sec) absolute value of the velocity (scale this per 1/3 second = 3fps)
    [angspeed,~] = getAngularSpeed3(centroid_smooth,1);
    TRACKS(i).angspeed = angspeed;
    
    %IDENTIFY HEAD AND TAIL, GET GRAYSCALE ASSOCIATED WITH THEM, THEN SMOOTH THEM, TOO
    [omega, head, tail, head_gs, tail_gs] = head_and_tail_kalman3(TRACKS(i));
    TRACKS(i).head = head;
    TRACKS(i).tail = tail;
    TRACKS(i).omega = omega;
    TRACKS(i).head_gs = head_gs;
    TRACKS(i).tail_gs = tail_gs;
    
    head_smooth = [movmean(TRACKS(i).head(:,1),3,'omitnan') movmean(TRACKS(i).head(:,2),3,'omitnan')];
    TRACKS(i).head_smooth = head_smooth;
    
    tail_smooth = [movmean(TRACKS(i).tail(:,1),3,'omitnan') movmean(TRACKS(i).tail(:,2),3,'omitnan')];
    TRACKS(i).tail_smooth = tail_smooth;
    
    %EXTRACT BOUTS OF FORWARD AND REVERSE MOVEMENT
    coherence_thresh = 90; %in degrees -- this requires that centroid vector and head or tail vector must be within this angle range of each other to be considered coherent motion
    speed_thresh = 0.02; %mm/sec -- this is the required speed to be considered moving
    [~, ~, ~, speed_smooth, forward, reverse] = getforwardreverse( centroid,head,tail,speed,coherence_thresh, speed_thresh);
    TRACKS(i).speed_smooth = speed_smooth;
    TRACKS(i).forward = forward;
    TRACKS(i).reverse = reverse;
    
    %DETERMINE WHEN ANIMAL WAS IN LAWN
    [headinlawn, centroidinlawn, tailinlawn, fullyinlawn] = countBlobsInOut(TRACKS(i), bg_struct);
    TRACKS(i).headinlawn = headinlawn;
    TRACKS(i).centroidinlawn = centroidinlawn;
    TRACKS(i).tailinlawn = tailinlawn;
    TRACKS(i).fullyinlawn = fullyinlawn;
    
    
    %ROAMING AND DWELLING
    %     trans = [0.995, 0.005; 0.07, 0.93]; %for HMM (Steve's numbers)
    %     emis = [0.96, 0.04; 0.07, 0.93];
    %     cutoff = 450;
    
    %my numbers from N2 OD2
    trans = [0.9645    0.0355; 0.0802    0.9198];
    emis =  [0.9790    0.0210; 0.5448    0.4552];
    cutoff = 35;
    x_offset = 2.5;
    binSize = 10*3;%10 seconds = 30 frames
    train = false;
    
    if TRACKS(i).age>=binSize %if the track length is longer than binSize, categories roaming and dwelling states
        [expSeq, expStates, ~, ~, actualratio] = getHMMStates3(TRACKS(i),binSize,TRACKS(i).centroidinlawn,cutoff,trans,emis,x_offset,train);
        TRACKS(i).roamdwell_2d = expSeq;
        TRACKS(i).roamdwell_hmm = expStates;
        TRACKS(i).actualratio = actualratio;
    else
        TRACKS(i).roamdwell_2d = NaN(TRACKS(i).age,1);
        TRACKS(i).roamdwell_hmm = NaN(TRACKS(i).age,1);
        TRACKS(i).actualratio = NaN(TRACKS(i).age,1);
    end
    
end

% REMOVE of dust or other objects that don't move
% USE MEAN SQUARED DISPLACEMENT
msds = {TRACKS.msd}';
avg_msd = cellfun(@nanmean,msds);
TRACKS = TRACKS(avg_msd>1e-5); %remove any tracks with inordinately low MSD

before_sec = 300; %5 minutes before
after_sec = 60; %1 minute after
%ENTER AND EXIT EVENTS, HEAD POKES
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

[TRACKS, POKE_INTS_GLOBAL, POKE_INTS_BY_TRACK, POKE_TRK_KEY, POKE_PEAK_IDX_GLOBAL, POKE_PEAK_IDX_BY_TRACK, POKE_DIST_MINSUBTRACT, POKE_IS_REV, POKE_RAD_DIST, POKE_EVHO_DIST, POKE_GS] = get_head_pokes2( TRACKS, bg_struct, INTS_OUT_BY_TRACK, OUT_INT_TRK_KEY );
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

startmin = 20; endmin = 60;
NUMWORMS = 1;
stat_int = [startmin*60*3 endmin*60*3]; %frames between 20 minutes and 60 minutes since the start of the assay
% stat_int = [1 7200]; %just for c4 and c5 on 05/08/18
curr_dir = pwd;
cd('..');
str = pwd ;
idx = strfind(str,'\') ;
foldername = str(idx(end)+1:end) ;
cd(curr_dir);
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
end

