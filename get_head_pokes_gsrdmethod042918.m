function [TRACKS, POKE_INTS_GLOBAL, POKE_INTS_BY_TRACK, POKE_TRK_KEY, POKE_PEAK_IDX_GLOBAL, POKE_PEAK_IDX_BY_TRACK, POKE_DIST_MINSUBTRACT, POKE_RAD_DIST, POKE_EVHO_DIST, POKE_GS] = get_head_pokes_gsrdmethod042918( TRACKS, bg_struct, INTS_OUT_BY_TRACK, OUT_INT_TRK_KEY )
%GET_HEAD_POKES_GSRDMETHOD042918.m This function identifies times that the worm pokes its
%head out of the lawn while the rest of the body remains in the lawn.
%This method looks specifically at peaks in the head grayscale value that
%occur at the lawn boundary.

%thresholds

% eh = bg_struct(TRACKS(1).bgvidindex(1)).ev_ho_crp_rel;
% [eh_cent_x, eh_cent_y, ~] = centroid(eh(:,1),eh(:,2));
% eh_rad_dist = diag(pdist2((repmat([eh_cent_x, eh_cent_y],length(eh),1)),eh)); %this is the distance of every point on the event horizon to the center of the lawn
% min_rd = min(eh_rad_dist)-10;
min_rd = 0; %is the ev_ho_distance stipulation enough perhaps?

gs_all = cell2mat({TRACKS(:).head_gs}');
hil_all = logical(cell2mat({TRACKS(:).headinlawn}'));
gsoutsidelawn = gs_all(~hil_all);
min_gs = nanmean(gsoutsidelawn)-3*nanstd(gsoutsidelawn); %minimum grayscale value to be considered outside the lawn is two standard deviations less than the mean grayscale value for the head outside the lawn.
% min_gs = 0.970; %maybe this just works better

min_prom = 5; %minimum peak prominence in the radial distance
time_tol = 5; %tolerance for peak overlap
dist_tol = 15; %distance tolerance = how close does the head have to be to the event horizon to be considered a head poke, 15 empirically seems to work

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
    smth_ed = movmean(EV_HO_DIST,5);
%     smth_ed = hampel(movmean(EV_HO_DIST,3),3);%moving average and remove jitters
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRAYSCALE / RADIAL-DISTANCE METHOD - look for peaks in radial
    % distance that occur close to peaks in head grayscale values (within a
    % certain distance to the event-horizon)
    
%     [~,gs_peaks, ~] = findpeaks(smth_gs,'MinPeakHeight',min_gs);
%     ehdist_at_gspeak = track.ev_ho_dist(gs_peaks);
%     gspeaksclosetoboundary = gs_peaks(ehdist_at_gspeak<dist_tol);
%     [closest_rd_peaks, rd_peak_intervals, rd_peak_height, idxlist_used] = findclosestpeakandinterval(gspeaksclosetoboundary,smth_rd,min_rd,min_prom,time_tol);
    
    %just try making it rd peaks where the grayscale is at least min_gs
    %that are within distance tolerance
    close_enough_to_boundary = smth_ed<dist_tol;
%     bright_enough_gs = smth_gs>=min_gs;
    placestocheckforrdpeaks = close_enough_to_boundary;% & bright_enough_gs;
    
    [peakheight,peak_centers,~,~, borders] = findpeaks_Elias(smth_rd,'MinPeakHeight',min_rd,'MinPeakProminence',min_prom,'WidthReference','halfheight'); %this modified version of the findpeaks method returns the borders of the peaks -- helpful for identifying intervals
    [closest_rd_peaks, rd_peak_intervals, rd_peak_height] = refine_peak_borders(smth_rd, peak_centers, borders, peakheight);
    [closest_rd_peaks,i_p] = setdiff(closest_rd_peaks,find(placestocheckforrdpeaks==0)); %select out only those rd peaks that occur in the accepted ranges by thresholds specified.
    rd_peak_intervals = rd_peak_intervals(i_p,:);
    rd_peak_height = rd_peak_height(i_p,:);
    
    left_flank = rd_peak_intervals(:,1); %beginning of head poke is the start of radial excursion outside the lawn
    right_flank = rd_peak_intervals(:,2);
    
    for idx = 1:length(closest_rd_peaks)
        potentialhp = closest_rd_peaks(idx);
        if isempty(left_flank(idx)) || isempty(right_flank(idx)) %don't include peaks at the beginning or end of a track, they can be unreliable
            continue;
        else
            %ADD A NEW POKE
            if ~logical(sum(potentialhp>banned_intervals(:,1) & potentialhp<banned_intervals(:,2))) && track.tailinlawn(potentialhp) %check that this head poke is not in one of the banned intervals (when the worm is OUT OF THE LAWN) AND that the tail is still in the lawn at the peak of the headpoke
                POKE_TRK_KEY = [POKE_TRK_KEY; trk];
                POKE_INTS_BY_TRACK = [POKE_INTS_BY_TRACK; left_flank(idx) right_flank(idx)];
                POKE_INTS_GLOBAL = [POKE_INTS_GLOBAL; left_flank(idx)+track.framesActive(1)-1 right_flank(idx)+track.framesActive(1)-1];
                HEAD_POKES(potentialhp)=1; %add this peak of poke to the HEAD_POKES vector
                POKE_PEAK_IDX_GLOBAL = [POKE_PEAK_IDX_GLOBAL; potentialhp+track.framesActive(1)-1];
                POKE_PEAK_IDX_BY_TRACK = [POKE_PEAK_IDX_BY_TRACK; potentialhp];
                POKE_DIST_MINSUBTRACT = [POKE_DIST_MINSUBTRACT; rd_peak_height(idx)-min_rd];
                POKE_RAD_DIST(struct_counter).trackid = trk;
                POKE_RAD_DIST(struct_counter).rad_dist = RAD_DIST(left_flank(idx):right_flank(idx))';
                POKE_EVHO_DIST(struct_counter).trackid = trk;
                POKE_EVHO_DIST(struct_counter).evho_dist = EV_HO_DIST(left_flank(idx):right_flank(idx))';
                POKE_GS(struct_counter).trackid = trk;
                POKE_GS(struct_counter).grayscale = HEAD_GS(left_flank(idx):right_flank(idx))';
                struct_counter = struct_counter+1;
            end
        end
    end
    
    TRACKS(trk).head_pokes = HEAD_POKES;
    TRACKS(trk).radial_dist = RAD_DIST; %save the non-smoothed version
    TRACKS(trk).ev_ho_dist = EV_HO_DIST;
end


end

%LOCAL FUNCTIONS
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