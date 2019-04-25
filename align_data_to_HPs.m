function ALIGNED_DATA = align_data_to_HPs( TRACKS, POKE_STRUCT, before_sec, after_sec, stat_int )
%ALIGN_DATA_TO_HPs.M This function takes in a set of tracks with pre-computed head poke events,
% and comes up with a set of intervals around head poke events.
% stacks the following data around these events:
% speed, angular speed, roaming and dwelling, head grayscale,
% ev_ho_dist, omegas, forward, reverse, fullyinlawn, centroidinlawn,
% headinlawn, tailinlawn <-- important to keep track of centroidinlawn in
% order to not average over data when the worm is outside of the lawn in
% case you want to know what worm behavior looks like before lawn-leaving.

after = (after_sec*3)-1; %in frames
before = (before_sec*3);
totallength = before+after+1;
ALIGNED_HP_IND = before+1;
ALL_HP = [POKE_STRUCT.POKE_TRK_KEY POKE_STRUCT.POKE_PEAK_IDX_BY_TRACK];
%only take those HPs that fall within the stat_int
ALL_HP = ALL_HP(POKE_STRUCT.POKE_PEAK_IDX_GLOBAL<=stat_int(2) & POKE_STRUCT.POKE_PEAK_IDX_GLOBAL>=stat_int(1),:);

ALL_VIDEONAME = cell(size(ALL_HP,1),totallength);
ALL_VIDEOFRAME = zeros(size(ALL_HP,1),totallength);
ALL_FRAMESACTIVE = zeros(size(ALL_HP,1),totallength);
ALL_BGVIDINDEX = zeros(size(ALL_HP,1),totallength);
ALL_CENTX = zeros(size(ALL_HP,1),totallength);
ALL_CENTY = zeros(size(ALL_HP,1),totallength);
ALL_HEADX = zeros(size(ALL_HP,1),totallength);
ALL_HEADY = zeros(size(ALL_HP,1),totallength);
ALL_TAILX = zeros(size(ALL_HP,1),totallength);
ALL_TAILY = zeros(size(ALL_HP,1),totallength);
ALL_SPEED = zeros(size(ALL_HP,1),totallength);
ALL_ANGSPEED = zeros(size(ALL_HP,1),totallength);
ALL_ROAMINGDWELLINGHMM = zeros(size(ALL_HP,1),totallength);
ALL_HEADGS = zeros(size(ALL_HP,1),totallength);
ALL_EVHODIST = zeros(size(ALL_HP,1),totallength);
ALL_RADDIST = zeros(size(ALL_HP,1),totallength);
ALL_HEADPOKES = zeros(size(ALL_HP,1),totallength);
ALL_HEADPOKE_REVERSAL = zeros(size(ALL_HP,1),totallength);
ALL_LAWNENTRIES = zeros(size(ALL_HP,1),totallength);
ALL_LAWNEXITS = zeros(size(ALL_HP,1),totallength);
ALL_FORWARD = zeros(size(ALL_HP,1),totallength);
ALL_REVERSE = zeros(size(ALL_HP,1),totallength);
ALL_OMEGA = zeros(size(ALL_HP,1),totallength);
ALL_CENTROIDINLAWN = zeros(size(ALL_HP,1),totallength);
ALL_HEADINLAWN = zeros(size(ALL_HP,1),totallength);
ALL_TAILINLAWN = zeros(size(ALL_HP,1),totallength);
ALL_FULLYINLAWN = zeros(size(ALL_HP,1),totallength);

HP_INTS = cell(size(ALL_HP,1),2); %hp index, which track the hp came from, which interval of that track the hp came from

tracks_fa = {TRACKS(:).framesActive}';
for i = 1:size(ALL_HP,1)
    track_ind = ALL_HP(i,1);                     % the track index to which this lawn-leaving event belongs
    track = TRACKS(track_ind);                      % the actual track
    hp_ind = ALL_HP(i,2);                      % the index of the leaving event (in track indices)
    
    curr_int_global = [track.framesActive(hp_ind)-before track.framesActive(hp_ind)+after];
    curr_int_inds_global = curr_int_global(1):curr_int_global(2);
    trk_key = [];                                                   % a list of the relevant tracks
    trk_ints = [];                                                  % a corresponding list of intervals in track indices
    trk_ints_align_idx = [];                                        % indices into curr_int_inds_global
    trk_key_per_hp = NaN(size(curr_int_inds_global));             % for every index -- which track the data came from
    trk_inds_per_hp = NaN(size(curr_int_inds_global));            % for every index -- which track inds they are
    for j = 1:size(tracks_fa,1) %get the track indices across all tracks that match the framesActive that surround a lawn-leaving event
        [~,i_ciig,i_tfa] = intersect(curr_int_inds_global,tracks_fa{j});
        if ~isempty(i_tfa)
            trk_key = [trk_key; j];
            trk_ints = [trk_ints; i_tfa(1) i_tfa(end)];
            trk_ints_align_idx = [trk_ints_align_idx; i_ciig(1) i_ciig(end)];
            trk_key_per_hp(i_ciig) = j;
            trk_inds_per_hp(i_ciig) = i_tfa;
        end
    end
    HP_INTS{i,1}=trk_key;
    HP_INTS{i,2}=trk_ints;
    
    %extract data from TRACKS using these indices and stack the data into ALIGNED_DATA
    vn = cell(1,totallength);
    vf = NaN(1,totallength);
    fa = NaN(1,totallength);
    bg = NaN(1,totallength);
    centx = NaN(1,totallength);
    centy = NaN(1,totallength);
    headx = NaN(1,totallength);
    heady = NaN(1,totallength);
    tailx = NaN(1,totallength);
    taily = NaN(1,totallength);
    speed = NaN(1,totallength);
    angspeed = NaN(1,totallength);
    roamdwellhmm = NaN(1,totallength);
    headgs = NaN(1,totallength);
    evhodist = NaN(1,totallength);
    raddist = NaN(1,totallength);
    headpokes = NaN(1,totallength);
    headpoke_reversal = NaN(1,totallength);
    lawnentries = NaN(1,totallength);
    lawnexits = NaN(1,totallength);
    forward = NaN(1,totallength);
    reverse = NaN(1,totallength);
    omega = NaN(1,totallength);
    centroidinlawn = NaN(1,totallength);
    headinlawn = NaN(1,totallength);
    tailinlawn = NaN(1,totallength);
    fullyinlawn = NaN(1,totallength);
    for k = 1:length(trk_key)
        vn(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))                 = TRACKS(trk_key(k)).videoname(trk_ints(k,1):trk_ints(k,2));
        vf(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))                 = TRACKS(trk_key(k)).videoframe(trk_ints(k,1):trk_ints(k,2));
        fa(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))                 = TRACKS(trk_key(k)).framesActive(trk_ints(k,1):trk_ints(k,2));
        bg(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))                 = TRACKS(trk_key(k)).bgvidindex(trk_ints(k,1):trk_ints(k,2));
        %         centx(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))            = TRACKS(trk_key(k)).centroid_smooth(trk_ints(k,1):trk_ints(k,2),1);
        %         centy(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))            = TRACKS(trk_key(k)).centroid_smooth(trk_ints(k,1):trk_ints(k,2),2);
        centx(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))              = TRACKS(trk_key(k)).centroid(trk_ints(k,1):trk_ints(k,2),1); %try the unsmoothed version
        centy(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))              = TRACKS(trk_key(k)).centroid(trk_ints(k,1):trk_ints(k,2),2);
        %         headx(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))            = TRACKS(trk_key(k)).head_smooth(trk_ints(k,1):trk_ints(k,2),1);
        %         heady(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))            = TRACKS(trk_key(k)).head_smooth(trk_ints(k,1):trk_ints(k,2),2);
        headx(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))              = TRACKS(trk_key(k)).head(trk_ints(k,1):trk_ints(k,2),1); %try the unsmoothed version
        heady(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))              = TRACKS(trk_key(k)).head(trk_ints(k,1):trk_ints(k,2),2);
        %         tailx(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))            = TRACKS(trk_key(k)).tail_smooth(trk_ints(k,1):trk_ints(k,2),1);
        %         taily(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))            = TRACKS(trk_key(k)).tail_smooth(trk_ints(k,1):trk_ints(k,2),2);
        tailx(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))              = TRACKS(trk_key(k)).tail(trk_ints(k,1):trk_ints(k,2),1); %try the unsmoothed version
        taily(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))              = TRACKS(trk_key(k)).tail(trk_ints(k,1):trk_ints(k,2),2);
        speed(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))              = TRACKS(trk_key(k)).speed_smooth(trk_ints(k,1):trk_ints(k,2));
        angspeed(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))           = TRACKS(trk_key(k)).angspeed(trk_ints(k,1):trk_ints(k,2));
        roamdwellhmm(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))       = TRACKS(trk_key(k)).roamdwell_hmm(trk_ints(k,1):trk_ints(k,2));
        headgs(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))             = TRACKS(trk_key(k)).head_gs(trk_ints(k,1):trk_ints(k,2));
        evhodist(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))           = TRACKS(trk_key(k)).ev_ho_dist(trk_ints(k,1):trk_ints(k,2));
        raddist(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))            = TRACKS(trk_key(k)).radial_dist(trk_ints(k,1):trk_ints(k,2));
        headpokes(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))          = TRACKS(trk_key(k)).head_pokes(trk_ints(k,1):trk_ints(k,2));
        
        %because i stupidly didn't make a logical array for tracks when
        %there are no hprev
        hpr = TRACKS(trk_key(k)).head_poke_reversal;
        if ~isempty(hpr)
            headpoke_reversal(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))  =  TRACKS(trk_key(k)).head_poke_reversal(trk_ints(k,1):trk_ints(k,2));
        end %otherwise, it remains NaNs
        
        lawnentries(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))        = TRACKS(trk_key(k)).lawn_entries(trk_ints(k,1):trk_ints(k,2));
        lawnexits(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))          = TRACKS(trk_key(k)).lawn_exits(trk_ints(k,1):trk_ints(k,2));
        forward(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))            = TRACKS(trk_key(k)).forward(trk_ints(k,1):trk_ints(k,2));
        reverse(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))            = TRACKS(trk_key(k)).reverse(trk_ints(k,1):trk_ints(k,2));
        omega(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))              = TRACKS(trk_key(k)).omega(trk_ints(k,1):trk_ints(k,2));
        centroidinlawn(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))     = TRACKS(trk_key(k)).centroidinlawn(trk_ints(k,1):trk_ints(k,2));
        headinlawn(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))         = TRACKS(trk_key(k)).headinlawn(trk_ints(k,1):trk_ints(k,2));
        tailinlawn(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))         = TRACKS(trk_key(k)).tailinlawn(trk_ints(k,1):trk_ints(k,2));
        fullyinlawn(trk_ints_align_idx(k,1):trk_ints_align_idx(k,2))        = TRACKS(trk_key(k)).fullyinlawn(trk_ints(k,1):trk_ints(k,2));
    end
    ALL_VIDEONAME(i,:) = vn;
    ALL_VIDEOFRAME(i,:) = vf;
    ALL_FRAMESACTIVE(i,:) = fa;
    ALL_BGVIDINDEX(i,:) = bg;
    ALL_CENTX(i,:) = centx;
    ALL_CENTY(i,:) = centy;
    ALL_HEADX(i,:) = headx;
    ALL_HEADY(i,:) = heady;
    ALL_TAILX(i,:) = tailx;
    ALL_TAILY(i,:) = taily;
    ALL_SPEED(i,:) = speed;
    ALL_ANGSPEED(i,:) = angspeed;
    ALL_ROAMINGDWELLINGHMM(i,:) = roamdwellhmm;
    ALL_HEADGS(i,:) = headgs;
    ALL_EVHODIST(i,:) = evhodist;
    ALL_RADDIST(i,:) = raddist;
    ALL_HEADPOKES(i,:) = headpokes;
    ALL_HEADPOKE_REVERSAL(i,:) = headpoke_reversal;
    ALL_LAWNENTRIES(i,:) = lawnentries;
    ALL_LAWNEXITS(i,:) = lawnexits;
    ALL_FORWARD(i,:) = forward;
    ALL_REVERSE(i,:) = reverse;
    ALL_OMEGA(i,:) = omega;
    ALL_CENTROIDINLAWN(i,:) = centroidinlawn;
    ALL_HEADINLAWN(i,:) = headinlawn;
    ALL_TAILINLAWN(i,:) = tailinlawn;
    ALL_FULLYINLAWN(i,:) = fullyinlawn;
end

ALIGNED_DATA = struct();
ALIGNED_DATA.HP_INTS = HP_INTS;
ALIGNED_DATA.ALIGNED_HP_IND = ALIGNED_HP_IND;
ALIGNED_DATA.VIDEONAME = ALL_VIDEONAME;
ALIGNED_DATA.VIDEOFRAME = ALL_VIDEOFRAME;
ALIGNED_DATA.FRAMESACTIVE = ALL_FRAMESACTIVE;
ALIGNED_DATA.BGVIDINDEX = ALL_BGVIDINDEX;
ALIGNED_DATA.CENTX = ALL_CENTX;
ALIGNED_DATA.CENTY = ALL_CENTY;
ALIGNED_DATA.HEADX = ALL_HEADX;
ALIGNED_DATA.HEADY = ALL_HEADY;
ALIGNED_DATA.TAILX = ALL_TAILX;
ALIGNED_DATA.TAILY = ALL_TAILY;
ALIGNED_DATA.SPEED = ALL_SPEED;
ALIGNED_DATA.ANGSPEED = ALL_ANGSPEED;
ALIGNED_DATA.ROAMINGDWELLINGHMM = ALL_ROAMINGDWELLINGHMM;
ALIGNED_DATA.HEADGS = ALL_HEADGS;
ALIGNED_DATA.EVHODIST = ALL_EVHODIST;
ALIGNED_DATA.RADDIST = ALL_RADDIST;
ALIGNED_DATA.HEADPOKES_ALL = ALL_HEADPOKES;
ALIGNED_DATA.HEADPOKES_REVERSAL = ALL_HEADPOKE_REVERSAL;
ALIGNED_DATA.LAWNENTRIES = ALL_LAWNENTRIES;
ALIGNED_DATA.LAWNEXITS = ALL_LAWNEXITS;
ALIGNED_DATA.FORWARD = ALL_FORWARD;
ALIGNED_DATA.REVERSE = ALL_REVERSE;
ALIGNED_DATA.OMEGA = ALL_OMEGA;
ALIGNED_DATA.CENTROIDINLAWN = ALL_CENTROIDINLAWN;
ALIGNED_DATA.HEADINLAWN = ALL_HEADINLAWN;
ALIGNED_DATA.TAILINLAWN = ALL_TAILINLAWN;
ALIGNED_DATA.FULLYINLAWN = ALL_FULLYINLAWN;

end

