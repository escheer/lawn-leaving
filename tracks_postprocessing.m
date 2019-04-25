function newAllTracks = tracks_postprocessing( allTracks, ev_ho_dict, bg_struct, pixpermm )
%TRACKS_POSTPROCESSING.m This function takes in the current tracks and the
%tracks that finished throughout the video and extracts a few more
%behavioral parameters from them.
%   Detailed explanation goes here

%% 
% remove jittery centroid positions to avoid any discontinuities
% calculate instantaneous speed from displacement
newAllTracks = struct(...
    'id', {}, ...
    'bgvidindex',{},...
    'x_shift',{},...
    'y_shift',{},...
    'centroid', {}, ...
    'bbox', {}, ...
    'cropworm',{}, ...
    'cropworm_orig',{},...
    'spline',{},...
    'curvature',{},...
    'posture_angle',{},...
    'worm',{},...
    'head',{},...
    'head_smooth',{},...
    'tail',{},...
    'tail_smooth',{},...
    'centroid_smooth',{},...
    'speed',{},...
    'speed_smooth',{},...
    'msd',{},...
    'msd_smooth',{},...
    'angspeed',{},...
    'roamdwell',{},...
    'actualratio',{},...
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
    %%% copy most fields directly to newAllTracks
    newAllTracks(i).id = allTracks(i).id;
    newAllTracks(i).x_shift = allTracks(i).x_shift;
    newAllTracks(i).y_shift = allTracks(i).y_shift;
    newAllTracks(i).centroid = allTracks(i).centroid;
    newAllTracks(i).bbox = allTracks(i).bbox;
    newAllTracks(i).cropworm = allTracks(i).cropworm;
    newAllTracks(i).cropworm_orig = allTracks(i).cropworm_orig;
    newAllTracks(i).spline = allTracks(i).spline;
    newAllTracks(i).worm = allTracks(i).worm;
    newAllTracks(i).curvature = allTracks(i).curvature;
    newAllTracks(i).posture_angle = allTracks(i).posture_angle;
    newAllTracks(i).age = allTracks(i).age;
    newAllTracks(i).consecutiveInvisibleCount = allTracks(i).consecutiveInvisibleCount;
    newAllTracks(i).totalVisibleCount = allTracks(i).totalVisibleCount;
    newAllTracks(i).framesActive = allTracks(i).framesActive;
    newAllTracks(i).videoname = allTracks(i).videoname;
    newAllTracks(i).videoframe = allTracks(i).videoframe;
    
    frames = allTracks(i).framesActive;
    bgvidindex = ev_ho_dict(frames);
    newAllTracks(i).bgvidindex = bgvidindex; %this can be used to index indto the bg_struct to access any of that information on a per frame basis
    %%%
    %SMOOTH CENTROID, GET MEAN SQUARED DISPLACEMENT
    centroid = newAllTracks(i).centroid;
    centroid_smooth = [movmean(centroid(:,1),3,'omitnan') movmean(centroid(:,2),3,'omitnan')];
    newAllTracks(i).centroid_smooth = centroid_smooth;
    
    [msd, msd_smth] = get_msd( centroid_smooth, pixpermm, 3 );
    newAllTracks(i).msd = msd;
    newAllTracks(i).msd_smooth = msd_smth;
    
    %GET SPEED AND ANGULAR SPEED
    set1 = centroid_smooth(1:end-3,:);
    set2 = centroid_smooth(4:end,:); %1 second ahead for calculation
    speed = [zeros(3,1); diag(pdist2(set1,set2))./pixpermm]; % D is the displacement between successive centroid positions in millimeters 
    newAllTracks(i).speed = speed; % (mm/sec) absolute value of the velocity (scale this per 1/3 second = 3fps)
    [angspeed,~] = getAngularSpeed2(centroid_smooth,3);
    newAllTracks(i).angspeed = angspeed;
    %IDENTIFY HEAD AND TAIL, GET GRAYSCALE ASSOCIATED WITH THEM, THEN SMOOTH THEM, TOO
    [omega, head, tail, head_gs, tail_gs] = head_and_tail_kalman(allTracks(i));
    newAllTracks(i).head = head;
    newAllTracks(i).tail = tail;
    newAllTracks(i).omega = omega;
    newAllTracks(i).head_gs = head_gs;
    newAllTracks(i).tail_gs = tail_gs;
    
    head_smooth = [movmean(newAllTracks(i).head(:,1),3,'omitnan') movmean(newAllTracks(i).head(:,2),3,'omitnan')];
    newAllTracks(i).head_smooth = head_smooth;
    
    tail_smooth = [movmean(newAllTracks(i).tail(:,1),3,'omitnan') movmean(newAllTracks(i).tail(:,2),3,'omitnan')];
    newAllTracks(i).tail_smooth = tail_smooth;
    
    %EXTRACT BOUTS OF FORWARD AND REVERSE MOVEMENT
    coherence_thresh = 90; %in degrees -- this requires that centroid vector and head or tail vector must be within this angle range of each other to be considered coherent motion
    speed_thresh = 0.02; %mm/sec -- this is the required speed to be considered moving
    [~, ~, ~, speed_smooth, forward, reverse] = getforwardreverse( centroid,head,tail,speed,coherence_thresh, speed_thresh);
    newAllTracks(i).speed_smooth = speed_smooth;
    newAllTracks(i).forward = forward;
    newAllTracks(i).reverse = reverse;
    
    %DETERMINE WHEN ANIMAL WAS IN LAWN
    [headinlawn, centroidinlawn, tailinlawn, fullyinlawn] = countBlobsInOut(newAllTracks(i), bg_struct);
    newAllTracks(i).headinlawn = headinlawn;
    newAllTracks(i).centroidinlawn = centroidinlawn;
    newAllTracks(i).tailinlawn = tailinlawn;
    newAllTracks(i).fullyinlawn = fullyinlawn;
    
    
    %ROAMING AND DWELLING
%     trans = [0.995, 0.005; 0.07, 0.93]; %for HMM (Steve's numbers)
%     emis = [0.96, 0.04; 0.07, 0.93];
%     cutoff = 450;

    trans = [0.967, 0.033; 0.043, 0.957];
    emis = [0.89, 0.11; 0.30, 0.70];
    cutoff = 1250;
    binSize = 5*3;%5 seconds = 15 frames
    
    [~, expStates, estTR, estE, actualratio] = getHMMStates(newAllTracks(i),binSize,centroidinlawn,cutoff,trans,emis);
    disp(estTR);
    disp(estE);
    newAllTracks(i).roamdwell = expStates;
    newAllTracks(i).actualratio = actualratio;
    
    %ENTER AND EXIT EVENTS
%     [enter_events, exit_events] = countCrossingTracks_after(newAllTracks(i), ev_ho_x, ev_ho_y);
%     newAllTracks(i).enter_events = enter_events;
%     newAllTracks(i).exit_events = exit_events;
end


end

