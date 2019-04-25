function TRACKS = tracks_postprocessing_JavierMods_032019( allTracks, ev_ho_dict, bg_struct, pixpermm, titlestr)
%TRACKS_POSTPROCESSING.m This function takes in the current tracks and the
%tracks that finished throughout the video and extracts a few more
%behavioral parameters from them.

%ROAMING AND DWELLING
%     trans = [0.995, 0.005; 0.07, 0.93]; %for HMM (Steve's numbers)
%     emis = [0.96, 0.04; 0.07, 0.93];
%     cutoff = 450;

%my numbers from N2 OD2 -- Javier will need his own!
trans = [0.9645    0.0355; 0.0802    0.9198];
emis =  [0.9790    0.0210; 0.5448    0.4552];
cutoff = 35;
x_offset = 2.5;
binSize = 10*3;%10 seconds = 30 frames
train = false;

TRACKS = struct(...
    'id', {}, ...
    'bgvidindex',{},...
    'x_shift',{},...
    'y_shift',{},...
    'centroid', {}, ...
    'bbox', {}, ...
    'cropworm',{}, ...
    'cropworm_orig',{},...
    'centroid_smooth',{},...
    'speed',{},...
    'msd',{},...
    'msd_smooth',{},...
    'angspeed',{},...
    'actualratio',{},...
    'roamdwell_hmm',{},...
    'roamdwell_2d',{},...
    'age', {}, ...
    'totalVisibleCount', {}, ...
    'consecutiveInvisibleCount', {}, ...
    'framesActive', {},...
    'videoname',{},...
    'videoframe',{},...
    'centroidinlawn',{});

outTrkIdx = 1;
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
    TRACKS(outTrkIdx).id = allTracks(i).id;
    TRACKS(outTrkIdx).x_shift = allTracks(i).x_shift(1:age);
    TRACKS(outTrkIdx).y_shift = allTracks(i).y_shift(1:age);
    TRACKS(outTrkIdx).centroid = allTracks(i).centroid(1:age,:);
    TRACKS(outTrkIdx).cropworm = allTracks(i).cropworm(1:age);
    TRACKS(outTrkIdx).cropworm_orig = allTracks(i).cropworm_orig(1:age);
    TRACKS(outTrkIdx).bbox = allTracks(i).bbox(1:age,:);
    TRACKS(outTrkIdx).age = age;
    TRACKS(outTrkIdx).consecutiveInvisibleCount = allTracks(i).consecutiveInvisibleCount;
    TRACKS(outTrkIdx).totalVisibleCount = age;
    TRACKS(outTrkIdx).framesActive = allTracks(i).framesActive(1:age);
    TRACKS(outTrkIdx).videoname = allTracks(i).videoname(1:age);
    TRACKS(outTrkIdx).videoframe = allTracks(i).videoframe(1:age);
    
    frames = TRACKS(outTrkIdx).framesActive;
    
    %to make this code extensible to previously tracked videos where
    %ev_ho_dict was not retained, check to see if the tracks already have
    %this field
    if ~isfield(allTracks(i),'bgvidindex')
        bgvidindex = ev_ho_dict(frames);
        TRACKS(outTrkIdx).bgvidindex = bgvidindex(1:age); %this can be used to index indto the bg_struct to access any of that information on a per frame basis
    else
        TRACKS(outTrkIdx).bgvidindex = allTracks(i).bgvidindex;
    end
    %%%
    %SMOOTH CENTROID, GET MEAN SQUARED DISPLACEMENT
    centroid = TRACKS(outTrkIdx).centroid;
    centroid_smooth = [movmean(centroid(:,1),3,'omitnan') movmean(centroid(:,2),3,'omitnan')];
    TRACKS(outTrkIdx).centroid_smooth = centroid_smooth;
    
    [msd, msd_smth] = get_msd( centroid_smooth, pixpermm, 3 );
    if nanmean(msd)<=1e-5 %if the msd is less than threshold, delete this track and move on. it is most likely dust, etc.
        TRACKS(outTrkIdx) = []; %dont increment outTrkIdx, so you re-use this index.
        continue;
    end
    TRACKS(outTrkIdx).msd = msd;
    TRACKS(outTrkIdx).msd_smooth = msd_smth;
    
    %GET SPEED AND ANGULAR SPEED
    set1 = centroid_smooth(1:end-3,:);
    set2 = centroid_smooth(4:end,:); %1 second ahead for calculation
    speed = [zeros(3,1); diag(pdist2(set1,set2))./pixpermm]; % D is the displacement between successive centroid positions in millimeters
    TRACKS(outTrkIdx).speed = speed;

	angspeed = getAngularSpeed_NavinMethod(centroid_smooth); %as of 2/27/19
    TRACKS(outTrkIdx).angspeed = angspeed;
    
    %DETERMINE WHEN ANIMAL CENTROID WAS IN LAWN
    track = TRACKS(outTrkIdx);
    centroidinlawn = zeros(track.age,1);
    for j = 1:track.age
        ev_ho = bg_struct(track.bgvidindex(j)).ev_ho_crp_rel; %get the event horizon that corresponds to this frame of the video
        x = ev_ho(:,1); y = ev_ho(:,2);
        centroidinlawn(j) = inpolygon(track.centroid(j,1),track.centroid(j,2),x,y);
    end
    TRACKS(outTrkIdx).centroidinlawn = centroidinlawn;

    if TRACKS(outTrkIdx).age>=binSize %if the track length is longer than binSize, categories roaming and dwelling states
        [expSeq, expStates, ~, ~, actualratio] = getHMMStates3(TRACKS(outTrkIdx),binSize,TRACKS(outTrkIdx).centroidinlawn,cutoff,trans,emis,x_offset,train);
        TRACKS(outTrkIdx).roamdwell_2d = expSeq;
        TRACKS(outTrkIdx).roamdwell_hmm = expStates;
        TRACKS(outTrkIdx).actualratio = actualratio;
    else
        TRACKS(outTrkIdx).roamdwell_2d = NaN(TRACKS(outTrkIdx).age,1);
        TRACKS(outTrkIdx).roamdwell_hmm = NaN(TRACKS(outTrkIdx).age,1);
        TRACKS(outTrkIdx).actualratio = NaN(TRACKS(outTrkIdx).age,1);
    end
    outTrkIdx = outTrkIdx+1;
end

end

