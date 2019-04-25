function [TRACKS, POKE_INTS_GLOBAL, POKE_INTS_BY_TRACK, POKE_TRK_KEY, POKE_PEAK_IDX_GLOBAL, POKE_PEAK_IDX_BY_TRACK, POKE_DIST_MINSUBTRACT, POKE_APPROACH_ANGLE, POKE_IS_FWD, POKE_IS_REV, POKE_SPEED, AVG_POKE_SPEED, POKE_RAD_DIST, POKE_EVHO_DIST, POKE_GS] = get_head_pokes3( TRACKS, bg_struct, INTS_OUT_BY_TRACK, OUT_INT_TRK_KEY, pixpermm )
%GET_HEAD_POKES3.m This function identifies times that the worm pokes its
%head out of the lawn while the rest of the body remains in the lawn.
%This method looks specifically at peaks in radial distance to define all
%headpokes and also demarcates a subcategory, headpoke-reversals, in which
%the headpoke is coupled to a backwards movement.

%thresholds
min_rd = 0;
min_prom = (5/112)*pixpermm; %minimum peak prominence in the radial distance (default = 5/112)
time_tol = 8; %tolerance for peak overlap (2.67 seconds)
dist_tol = (15/112)*pixpermm; %distance tolerance = how close does the head have to be to the event horizon to be considered a head poke, 15 empirically seems to work
speed_thresh = 0.02;

POKE_INTS_GLOBAL = []; %keeps track of the head poke intervals in the whole video, indexed by global frame numbers
POKE_INTS_BY_TRACK = []; %same but indices are from start of track = 1
POKE_TRK_KEY = []; %corresponding to the last two vectors -- which track do they refer to among TRACKS?
POKE_PEAK_IDX_GLOBAL = [];
POKE_PEAK_IDX_BY_TRACK = [];
POKE_DIST_MINSUBTRACT = [];
POKE_IS_REV = [];
POKE_IS_FWD = [];
POKE_APPROACH_ANGLE = [];
POKE_SPEED = struct();
AVG_POKE_SPEED = [];
POKE_RAD_DIST = struct();
POKE_EVHO_DIST = struct();
POKE_GS = struct();
struct_counter = 1;

for trk = 1:length(TRACKS)
    %     disp(trk);
    donewithtrack = false;
    noreversals = false;
    
    track = TRACKS(trk);
    banned_intervals = INTS_OUT_BY_TRACK(OUT_INT_TRK_KEY==trk,:);
    if isempty(banned_intervals)
        banned_intervals = [-10 -10]; %an interval that will never be a problem
    end
    
    EV_HO_DIST = NaN(track.age,1);
    RAD_DIST = NaN(track.age,1);
    HEAD_POKES = false(track.age,1);
    HEAD_POKE_REV = false(track.age,1);
    HEAD_POKE_ANGLE = NaN(track.age,1); %NaNs everywhere, and angle of approach only at head pokes
    
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
    RAD_DIST = fillmissing(RAD_DIST,'linear');
    EV_HO_DIST = fillmissing(EV_HO_DIST,'linear');
    
    %smooth data
    smth_rd = RAD_DIST;
    smth_ed = movmean(EV_HO_DIST,5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RADIAL-DISTANCE METHOD - look for peaks in radial
    % distance that occur close to peaks in head grayscale values (within a
    % certain distance to the event-horizon) when the animal is moving
    % forward;
    
    % conditions to meet to be considered at head poke
    close_enough_to_boundary = smth_ed<dist_tol;
%     moving_forward = track.forward;
    conditions_met = close_enough_to_boundary;% & moving_forward;
    
    [peakheight,peak_centers,~,~, borders] = findpeaks_Elias(smth_rd,'MinPeakHeight',min_rd,'MinPeakProminence',min_prom,'WidthReference','halfheight'); %this modified version of the findpeaks method returns the borders of the peaks -- helpful for identifying intervals
    [rd_peaks_ALL, rd_peak_intervals_ALL, rd_peak_height_ALL] = refine_peak_borders(smth_rd, peak_centers, borders, peakheight);
    [rd_peaks_ALL,i_p] = setdiff(rd_peaks_ALL,find(conditions_met==0)); %select out only those rd peaks that occur in the accepted ranges by thresholds specified.
    rd_peak_intervals_ALL = rd_peak_intervals_ALL(i_p,:);
    rd_peak_height_ALL = rd_peak_height_ALL(i_p);
   
    idx_to_remove = [];
    % make sure that the worm is moving forward before the head poke
    for k = 1:size(rd_peak_intervals_ALL,1)
        int_idx = rd_peak_intervals_ALL(k,1):rd_peaks_ALL(k);
        speed_during_int = track.speed(int_idx);
        if sum(abs(speed_during_int)<speed_thresh)/length(speed_during_int)>0.50 %if the animal is not moving for less than 50% of interval preceding HP, delete it
            idx_to_remove = [idx_to_remove; rd_peaks_ALL(k)];
        end
    end
    
    % if there are multiple radial distance peaks close together while the
    % head is outside the lawn, just choose the biggest one.
    head_out_ints = get_intervals( ~track.headinlawn, 1 );
    for k = 1:size(head_out_ints,1)
        rdpeaks_headout_tf = ismember(rd_peaks_ALL,head_out_ints(k,1):head_out_ints(k,2)); %1 or 0 for every peak
        if sum(rdpeaks_headout_tf)>1
            peakheadoutidx = find(rdpeaks_headout_tf); %indices of peaks outside the lawn in a contiguous interval (can index into rd_peaks_ALL)
            rdpeaks_headout_height = rd_peak_height_ALL(peakheadoutidx); %height of said peaks
            [~,tmp_idx] = max(rdpeaks_headout_height);
            idx_to_remove = [idx_to_remove ; rd_peaks_ALL(peakheadoutidx(setdiff(1:length(rdpeaks_headout_height),tmp_idx)))]; %these are the peaks to get rid of (not the heighest one)
        end
    end
    % remove any duplicates
    idx_to_remove = unique(idx_to_remove);
    %find these indices in the peaks lists
    to_remove_peak_idx = ismember(rd_peaks_ALL,idx_to_remove);
    rd_peaks_ALL(to_remove_peak_idx) = []; 
    rd_peak_intervals_ALL(to_remove_peak_idx,:) = [];
    rd_peak_height_ALL(to_remove_peak_idx) = [];

    left_flank_ALL = rd_peak_intervals_ALL(:,1); %beginning of head poke is the start of radial excursion outside the lawn
    right_flank_ALL = rd_peak_intervals_ALL(:,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REVERSAL METHOD - look for reversals that occur close to the
    % lawn boundary -- check if there was a peak in head grayscale, too
    reversalints = get_intervals( track.reverse, 1 ); %get intervals when worm was reversing
    
    if ~isempty(reversalints)
        reversalstarts = reversalints(:,1);
    else
        reversalstarts = [];
        noreversals = true;
    end
    if isempty(reversalstarts) || isempty(rd_peaks_ALL) %if there are no reversals close enough to the boundary for this track, continue to the next track.
        donewithtrack = true;
    end
    
    %CATEGORIZE HEAD POKES
    if ~donewithtrack %if there ARE reversals close enough to the boundary
        if noreversals
            idx_select = [];
        else
            [~, idx_select] = findclosestreversaltordpeak(reversalstarts,rd_peaks_ALL,time_tol);
        end

        for idx = 1:length(rd_peaks_ALL)
            potentialhp = rd_peaks_ALL(idx);
            if isempty(left_flank_ALL(idx)) || isempty(right_flank_ALL(idx)) %don't include peaks at the beginning or end of a track, they can be unreliable
                continue;
            else
                %ADD A NEW POKE
                if ~logical(sum(potentialhp>banned_intervals(:,1) & potentialhp<banned_intervals(:,2))) && track.tailinlawn(potentialhp) %check that this head poke is not in one of the banned intervals (when the worm is OUT OF THE LAWN) AND that the tail is still in the lawn at the peak of the headpoke
                    POKE_TRK_KEY = [POKE_TRK_KEY; trk];
                    POKE_INTS_BY_TRACK = [POKE_INTS_BY_TRACK; left_flank_ALL(idx) right_flank_ALL(idx)];
                    POKE_INTS_GLOBAL = [POKE_INTS_GLOBAL; left_flank_ALL(idx)+track.framesActive(1)-1 right_flank_ALL(idx)+track.framesActive(1)-1];
                    HEAD_POKES(potentialhp)=1; %add this peak of poke to the HEAD_POKES vector
                    POKE_PEAK_IDX_GLOBAL = [POKE_PEAK_IDX_GLOBAL; potentialhp+track.framesActive(1)-1];
                    POKE_PEAK_IDX_BY_TRACK = [POKE_PEAK_IDX_BY_TRACK; potentialhp];
                    POKE_DIST_MINSUBTRACT = [POKE_DIST_MINSUBTRACT; rd_peak_height_ALL(idx)-min_rd];
                    POKE_IS_FWD = [POKE_IS_FWD; true];
                    POKE_IS_REV = [POKE_IS_REV; false];
                    %get poke approach angle
                    eh = bg_struct(track.bgvidindex(potentialhp)).ev_ho_crp_rel; %get event horizon for this head poke
                    angle = getpokeapproachangle(track.head,track.centroid,eh,[left_flank_ALL(idx) potentialhp]);
                    POKE_APPROACH_ANGLE = [POKE_APPROACH_ANGLE ; angle];
                    HEAD_POKE_ANGLE(potentialhp) = angle;
                    AVG_POKE_SPEED = [AVG_POKE_SPEED; nanmean(track.speed(left_flank_ALL(idx):potentialhp))]; %average speed of the radial excursion (until peak)
                    POKE_SPEED(struct_counter).trackid = trk;
                    POKE_SPEED(struct_counter).speed = track.speed(left_flank_ALL(idx):right_flank_ALL(idx))';
                    POKE_RAD_DIST(struct_counter).trackid = trk;
                    POKE_RAD_DIST(struct_counter).rad_dist = RAD_DIST(left_flank_ALL(idx):right_flank_ALL(idx))';
                    POKE_EVHO_DIST(struct_counter).trackid = trk;
                    POKE_EVHO_DIST(struct_counter).evho_dist = EV_HO_DIST(left_flank_ALL(idx):right_flank_ALL(idx))';
                    POKE_GS(struct_counter).trackid = trk;
                    POKE_GS(struct_counter).grayscale = HEAD_GS(left_flank_ALL(idx):right_flank_ALL(idx))';
                    struct_counter = struct_counter+1;
                    if ismember(idx,idx_select) %this is also a HP_REVERSAL
                        HEAD_POKE_REV(potentialhp)=true;
                        POKE_IS_REV(end)=true; %change this to a true
                        POKE_IS_FWD(end)=false;
                    end
                end
            end
        end
    end
    TRACKS(trk).head_pokes = HEAD_POKES;
    TRACKS(trk).head_poke_forward = HEAD_POKES & ~HEAD_POKE_REV; %new field for just the forward head pokes
    TRACKS(trk).head_poke_reversal = HEAD_POKE_REV;
    TRACKS(trk).radial_dist = RAD_DIST; %save the non-smoothed version
    TRACKS(trk).ev_ho_dist = EV_HO_DIST;
    TRACKS(trk).head_poke_angle = HEAD_POKE_ANGLE;
end
POKE_IS_FWD = logical(POKE_IS_FWD);
POKE_IS_REV = logical(POKE_IS_REV);
end

%LOCAL FUNCTIONS
function [rd_peak_select, idx_select] = findclosestreversaltordpeak(rev_idx,rd_peaks,timetolerance)
dist_mat = pdist2(rd_peaks,rev_idx); %time separation of rev_idx and peaks in rd
[min_dist , idx] = min( dist_mat,[],1 ); %get the closest peaks
select = min_dist < timetolerance;  %select out those close enough together based on the tolerance we specified above.
% idxlist_used = find(select); %which of the original reversals have a peak in radial distance close enough
idx_select = unique(idx(select)); %the same peak may be chosen multiple times -- ensure that it is listed only once.
rd_peak_select = rd_peaks(idx_select); %get the actual matches (in rd)
end

function angle = getpokeapproachangle(head,cent,eh,poke_int)
y = eh(:,2);
x = eh(:,1);

poke_head = head(poke_int(1):poke_int(end),:);

% 1a. look for intersections of the head path with the event horizon. if
% one exists, add it to the event horizon points
head_crossings = InterX([x y]',poke_head')';
if ~isempty(head_crossings) && size(head_crossings,1)==1 && ~ismember(head_crossings,[x y],'rows')
    x = [x; head_crossings(1)]; y = [y; head_crossings(2)]; %add the intersection point to event horizon
    [x,y] = sortPointsCw(x,y); %make sure the points are in order (clockwise)
    closest_eh_point = head_crossings;
    [~,which_eh_point] = min(sqrt((x-closest_eh_point(1)).^2 + (y-closest_eh_point(2)).^2)); %find out the index of the head crossing in the new event horizon
else
    [eh_points,head_dist,~] = distance2curve([x y],poke_head);
    % 1b. look for the closest point in the head trajectory leading up to HP to
    % a point on the event horizon (if there was a head intersection, this
    % should be that point)
    [~,h_i] = min(head_dist);
    closest_eh_point = eh_points(h_i,:);
    [~,which_eh_point] = min(sqrt((x-closest_eh_point(1)).^2 + (y-closest_eh_point(2)).^2));
end
slope = gradient(y,x);

% 2. what is the tangent line to that point?
x_range = (x(which_eh_point)-5:1:x(which_eh_point)+5)';
tangent = (x_range-x(which_eh_point))*slope(which_eh_point)+y(which_eh_point);
%find a local segment of this tangent that does not loop over the range of x 
t_vec = [x_range(end) tangent(end)]-[x_range(1) tangent(1)]; t_vec = [t_vec 0];

% 3. find the centroid vector through the points in poke_int
smth_window = 3;
c_smth = [movmean(cent(:,1),smth_window,'omitnan') movmean(cent(:,2),smth_window,'omitnan')];
poke_cent = c_smth(poke_int(1):poke_int(end),:);
% c_vec = poke_cent(end,:)-poke_cent(1,:); c_vec = [c_vec 0]; %centroid vector formulated as the last centroid position minus the first
c_vec = mean(diff(poke_cent,1,1),1); c_vec = [c_vec 0]; %alternate formulation - average centroid during lead up to head poke

% 4. calculate the angle between the centroid vector and the closest point
% tangent line
angle = abs(atan2d(norm(cross(c_vec,t_vec)),dot(c_vec,t_vec)));
if angle>90 %may need to subtract angle from 180 to get the angle between 0 and 90
    angle = 180-angle;
end

% disp(['Angle of Approach is ... ' num2str(angle)]);
% h=figure(); hold on;
% plot(x,y); axis equal;
% plot(poke_head(:,1),poke_head(:,2));
% plot(poke_cent(:,1),poke_cent(:,2));
% scatter(poke_cent(1,1),poke_cent(1,2),1,'b+');
% scatter(poke_cent(end,1),poke_cent(end,2),1,'r+');
% scatter(closest_eh_point(1),closest_eh_point(2),2,'g+');
% plot(x_range,tangent);
% pause(); close(h);
end


function [closestpeaks, peakintervals, peakheight, idxlist_used] = findclosestpeakandinterval(idxlist,timevec,minheight,minprom,timetolerance) %local function that finds the closest peak (and their intervals) in timevec to each index in idxlist

[peak_height,locs,~,~, locs_borders] = findpeaks_Elias(timevec,'MinPeakHeight',minheight,'MinPeakProminence',minprom,'WidthReference','halfheight'); %this modified version of the findpeaks method returns the borders of the peaks -- helpful for identifying intervals
dist_mat = pdist2(locs,idxlist); %time separation of peaks idxlist and peaks in timevec
[min_dist , idx] = min( dist_mat,[],1 ); %get the closest peaks
select = min_dist < timetolerance;  %select out those close enough together based on the tolerance we specified above.
idxlist_used = find(select); %which of the original reversals have a peak in radial distance close enough
idx_select = unique(idx(select)); %the same peak may be chosen multiple times -- ensure that it is listed only once.

unique_matches = locs(idx_select); %get the actual matches (in rd)
borders = locs_borders(idx_select,:);
peakheight = peak_height(idx_select);

[closestpeaks, peakintervals, peakheight] = refine_peak_borders(timevec, unique_matches, borders, peakheight);


end

function [peaks, peak_intervals, peak_height] = refine_peak_borders(timevec, peak_centers, borders, peakheight)
%refine these boundaries by looking for intersections of the line
%of half-maximum height for each peak with the peak on either side. (by
%interpolation).
peaks = zeros(length(peak_centers),1);
peak_intervals = zeros(length(peak_centers),2);
peak_height = zeros(length(peak_centers),1);

halfheight = peakheight./2;
for i = 1:length(halfheight)
    peak_center = peak_centers(i);
    leftborder = borders(i,1); rightborder = borders(i,2);
    halfheight_line = [leftborder rightborder ; halfheight(i) halfheight(i)]; %first row is x values of borders, second row is the half height (horizontal line)
    chunk = [leftborder:rightborder; timevec(leftborder:rightborder)']; %chunk of the smth vector in which to look for crossings
    crossings = InterX(halfheight_line,chunk); %find the intersections
    left_flank = leftborder; right_flank = rightborder; %initialize
    
    if ~isempty(crossings) %if there are intersections within the borders, use those instead, since they are more refined.
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
    peaks(i) = peak_center;
    peak_intervals(i,:) = [left_flank right_flank];
    peak_height(i) = peakheight(i);
end

end

