function BORDER_DATA = extract_intervals_at_boundary( TRACKS, stat_int, ehdist_thresh, min_duration )
%EXTRACT_INTERVALS_AT_BOUNDARY.M This function takes in a set of tracks and
%extracts intervals when the animal is within a threshold distance of the
%food boundary for at least min_duration frames. 

BORDER_DATA = struct();
counter = 1;
for i = 1:length(TRACKS)
    track = TRACKS(i);
    atborder = track.ev_ho_dist < ehdist_thresh; %get stretches where the animal is within threshold distance of the border
    atborder_ints = get_intervals(atborder,1)+track.framesActive(1)-1; %get the start and stop indices of these intervals (put temporarily into global indices)
    if isempty(atborder_ints)
        continue;
    end
    atborder_ints = atborder_ints(atborder_ints(:,1)>=stat_int(1) & atborder_ints(:,2)>=stat_int(1) & atborder_ints(:,1)<=stat_int(2) & atborder_ints(:,2)<=stat_int(2),:); %only consider intervals in the stat_int
    atborder_ints = atborder_ints-track.framesActive(1)+1; %go back to track indices
    atborder_dur = atborder_ints(:,2)-atborder_ints(:,1);
    atborder_ints = atborder_ints(atborder_dur>min_duration,:); %only consider stretches lasting at least min_duration
        
    for k = 1:size(atborder_ints,1) %go through each interval at the border and make a separate track out of it to add to the border data.
        idx = atborder_ints(k,1):atborder_ints(k,2);
        for f = fieldnames(track)'
%             f
            if strcmp(f{1},'end1') || strcmp(f{1},'end2') || strcmp(f{1},'end1_g') || strcmp(f{1},'end2_g')
                continue;
            end
            if length(track.(f{1}))==1
                BORDER_DATA(counter).(f{1}) = track.(f{1});
            else
                if isempty(track.(f{1}))
                    BORDER_DATA(counter).(f{1}) = [];
                else
                    BORDER_DATA(counter).(f{1}) = track.(f{1})(idx);
                end
            end
        end
        counter = counter+1;
    end
end


end

