function [ SPEED, SPEED_ON_FOOD, MEAN_SPEED_ON_FOOD, SPEED_OFF_FOOD, MEAN_SPEED_OFF_FOOD, ANGSPEED, ROAMDWELL_2D, ROAMDWELL_HMM, FRAC_ROAMING_HMM ] = get_SPEED_ANGSPEED_RD_during_interval( TRACKS, stat_int, OK_frames_inlawn )
%get_SPEED_ANGSPEED_RD_during_interval.m This function collates speed,
%angular speed, and roaming and dwelling information from all tracks and
%selects out the data that falls into the acceptable interval.

SPEED = NaN(7200,1);
SPEED_ON_FOOD = NaN(7200,1);
SPEED_OFF_FOOD = NaN(7200,1);
ANGSPEED = NaN(7200,1);
ROAMDWELL_2D = NaN(7200,1);
ROAMDWELL_HMM = NaN(7200,1);

for trk = 1:length(TRACKS)
    track = TRACKS(trk);
    [OK_trk_idx, OK_int_idx] = ismember(track.framesActive,stat_int(1)+1:stat_int(2)); %have to add +1 to make it 7200 long
    OK_trk_idx = find(OK_trk_idx);
    OK_int_idx = OK_int_idx(OK_int_idx~=0);
    if ~isempty(OK_trk_idx)
        SPEED(OK_int_idx) = track.speed(OK_trk_idx);
        ANGSPEED(OK_int_idx) = track.angspeed(OK_trk_idx);
        ROAMDWELL_2D(OK_int_idx) = track.roamdwell_2d(OK_trk_idx);
        ROAMDWELL_HMM(OK_int_idx) = track.roamdwell_hmm(OK_trk_idx);
        [OK_trk_inlawn_idx, OK_int_inlawn_idx] = ismember(track.framesActive,OK_frames_inlawn); %get the speed when the animal is only ON FOOD
        OK_trk_inlawn_idx = find(OK_trk_inlawn_idx);
        OK_int_inlawn_idx = OK_int_inlawn_idx(OK_int_inlawn_idx~=0);
        if ~isempty(OK_trk_inlawn_idx)
            SPEED_ON_FOOD(OK_int_inlawn_idx) = track.speed(OK_trk_inlawn_idx);
        end
        OK_trk_outlawn_idx = setdiff(OK_trk_idx,OK_trk_inlawn_idx);
        OK_int_outlawn_idx = setdiff(OK_int_idx,OK_int_inlawn_idx);
        if ~isempty(OK_trk_outlawn_idx)
            SPEED_OFF_FOOD(OK_int_outlawn_idx) = track.speed(OK_trk_outlawn_idx);
        end
    else
        continue;
    end
end
MEAN_SPEED_ON_FOOD = nanmean(SPEED_ON_FOOD);
MEAN_SPEED_OFF_FOOD = nanmean(SPEED_OFF_FOOD);
FRAC_ROAMING_HMM = nansum(ROAMDWELL_HMM==2)/sum(~isnan(ROAMDWELL_HMM));

end

