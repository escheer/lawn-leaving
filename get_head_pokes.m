function [TRACKS, POKE_INTS_GLOBAL, POKE_INTS_BY_TRACK, POKE_TRK_KEY, POKE_PEAK_IDX_GLOBAL, POKE_PEAK_IDX_BY_TRACK, POKE_DIST_MINSUBTRACT, POKE_RAD_DIST, POKE_EVHO_DIST, POKE_GS] = get_head_pokes( TRACKS, bg_struct, min_gs, INTS_OUT_BY_TRACK, OUT_INT_TRK_KEY )
%GET_HEAD_POKES.m This function identifies times that the worm pokes its
%head out of the lawn while the rest of the body remains in the lawn.

eh = bg_struct(TRACKS(1).bgvidindex(1)).ev_ho_crp_rel;
[eh_cent_x, eh_cent_y, ~] = centroid(eh(:,1),eh(:,2));
eh_rad_dist = diag(pdist2((repmat([eh_cent_x, eh_cent_y],length(eh),1)),eh)); %this is the distance of every point on the event horizon to the center of the lawn
min_rd = min(eh_rad_dist)-10;
tol = 3; %tolerance for peak overlap

POKE_INTS_GLOBAL = []; %keeps track of the head poke intervals in the whole video, indexed by global frame numbers
POKE_INTS_BY_TRACK = []; %same but indices are from start of track = 1
POKE_TRK_KEY = []; %corresponding to the last two vectors -- which track do they refer to among TRACKS?
POKE_PEAK_IDX_GLOBAL = [];
POKE_PEAK_IDX_BY_TRACK = [];
POKE_DIST_MINSUBTRACT = [];
POKE_RAD_DIST = struct();
POKE_EVHO_DIST = struct();
POKE_GS = struct();
struct_counter = 1;

for trk = 1:length(TRACKS)
    disp(trk);
    track = TRACKS(trk);
    banned_intervals = INTS_OUT_BY_TRACK(OUT_INT_TRK_KEY==trk,:);
    if isempty(banned_intervals)
        banned_intervals = [-10 -10]; %an interval that will never be a problem
    end
%         if trk == 7
%             disp('debug');
%         end
    EV_HO_DIST = zeros(track.age,1);
    RAD_DIST = zeros(track.age,1);
    HEAD_POKES = zeros(track.age,1);
    for j = 1:track.age
        eh = bg_struct(track.bgvidindex(j)).ev_ho_crp_rel;
        [eh_cent_x, eh_cent_y, ~] = centroid(eh(:,1),eh(:,2));
        head = track.head(j,:);
        RAD_DIST(j) = pdist2(head,[eh_cent_x eh_cent_y]); %find head distance from the center of the lawn
        [~,dist,~] = distance2curve(eh,head); %find head distance to the nearest point on event horizon
        EV_HO_DIST(j) = dist;
    end
    tmp = EV_HO_DIST(~track.headinlawn);%make distance negative if its outside the event horizon
    EV_HO_DIST(~track.headinlawn) = -1*tmp;
    HEAD_GS = track.head_gs;
    %fill in missing data, this improves analysis
    HEAD_GS = fillmissing(HEAD_GS,'linear');
    RAD_DIST = fillmissing(RAD_DIST,'linear');
    EV_HO_DIST = fillmissing(EV_HO_DIST,'linear');
    %smooth data
    smth_gs = movmean(HEAD_GS,3);
    smth_rd = movmean(RAD_DIST,5);
    smth_ed = hampel(movmean(EV_HO_DIST,5),3);%moving average and remove jitters
    
    %now find peaks of these two quantities and impose the thresholds
    %established above
    [~,gs_locs, ~] = findpeaks(smth_gs,'MinPeakHeight',min_gs);
    [rd_peak_height,rd_locs,~,~, rd_locs_borders] = findpeaks_Elias(smth_rd,'MinPeakHeight',min_rd,'MinPeakProminence',5,'WidthReference','halfheight');
    dist_mat = pdist2(rd_locs,gs_locs);
    [min_dist_rd , idx_rd ] = min( dist_mat,[],1 );
    rd_select = min_dist_rd < tol;
    ofinterest = unique(idx_rd(rd_select));
    unique_matches_in_rd = rd_locs(ofinterest);
    corresponding_borders = rd_locs_borders(ofinterest,:);
    corresponding_peak_height = rd_peak_height(ofinterest);
    
    %now refine these boundaries by looking for intersections of the line
    %of half-maximum height for each peak with the peak on either side. (by
    %interpolation).
    halfheight = corresponding_peak_height./2;
    for i = 1:length(halfheight)
        peak_center = unique_matches_in_rd(i);
        leftborder = corresponding_borders(i,1); rightborder = corresponding_borders(i,2);
        %         leftborder = 1; rightborder = track.age;
        halfheight_line = [leftborder rightborder ; halfheight(i) halfheight(i)]; %first row is x values of borders, second row is the half height (horizontal line)
        rd_chunk = [leftborder:rightborder; smth_rd(leftborder:rightborder)']; %chunk of the smth_rd vector in which to look for crossings
        crossings = InterX(halfheight_line,rd_chunk);
        left_flank = leftborder; right_flank = rightborder;
        
        if ~isempty(crossings) %if there are intersections within the borders, use those instead
            x_vals = crossings(1,:);
            lessthan = x_vals(x_vals<peak_center);
            if ~isempty(lessthan)
                test = round(max(lessthan)); %greatest member less than peak center
                if test > left_flank
                    left_flank = test;
                end
            end
            greaterthan = x_vals(x_vals>peak_center);
            if ~isempty(greaterthan)
                test = round(min(greaterthan)); %least member greater than peak center
                if test < right_flank
                    right_flank = test;
                end
            end
        end
        
        if isempty(left_flank) || isempty(right_flank) %don't include peaks at the beginning or end of a track, they can be unreliable
            continue;
        else
            %ADD A NEW POKE
            if ~logical(sum(peak_center>banned_intervals(:,1) & peak_center<banned_intervals(:,2))) && track.tailinlawn(peak_center) %check that this head poke is not in one of the banned intervals (when the worm is OUT OF THE LAWN) AND that the tail is still in the lawn at the peak of the headpoke
                POKE_TRK_KEY = [POKE_TRK_KEY; trk];
                POKE_INTS_BY_TRACK = [POKE_INTS_BY_TRACK; left_flank right_flank];
                POKE_INTS_GLOBAL = [POKE_INTS_GLOBAL; left_flank+track.framesActive(1)-1 right_flank+track.framesActive(1)-1];
                HEAD_POKES(peak_center)=1; %add this peak of poke to the HEAD_POKES vector
                POKE_PEAK_IDX_GLOBAL = [POKE_PEAK_IDX_GLOBAL; peak_center+track.framesActive(1)-1];
                POKE_PEAK_IDX_BY_TRACK = [POKE_PEAK_IDX_BY_TRACK; peak_center];
                POKE_DIST_MINSUBTRACT = [POKE_DIST_MINSUBTRACT; corresponding_peak_height(i)-min_rd];
                POKE_RAD_DIST(struct_counter).trackid = trk;
                POKE_RAD_DIST(struct_counter).rad_dist = smth_rd(left_flank:right_flank)';
                POKE_EVHO_DIST(struct_counter).trackid = trk;
                POKE_EVHO_DIST(struct_counter).evho_dist = smth_ed(left_flank:right_flank)';
                POKE_GS(struct_counter).trackid = trk;
                POKE_GS(struct_counter).grayscale = smth_gs(left_flank:right_flank)';
                struct_counter = struct_counter+1;
            end
        end
    end
   
    TRACKS(trk).head_pokes = HEAD_POKES;
    TRACKS(trk).radial_dist = RAD_DIST; %save the non-smoothed version
    TRACKS(trk).ev_ho_dist = EV_HO_DIST;
end


end

