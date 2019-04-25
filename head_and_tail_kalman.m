function [omega, head, tail, head_gs, tail_gs] = head_and_tail_kalman( WORMTRACK )
%head_and_tail_kalman.m This function tracks each endpoint of the spline
%independently with kalman filters.

origPlayer = vision.VideoPlayer('Position', [120, 160, 100, 100]);

lastgood = 1;
justbroke = 0;
badcount = 0;

[tracks,tracksThatLeft] =  initializeTracks();
nextId = 1;

centroid_dist = NaN(WORMTRACK.totalVisibleCount,1);
possible_fixes = [];
for j = 1:WORMTRACK.totalVisibleCount
    disp(['idx: ' num2str(j)]);

    
    x_off = double(WORMTRACK.bbox(j,1))-5;
    y_off = double(WORMTRACK.bbox(j,2))-5;
    
    EvWorm = WORMTRACK.worm{j};
    
    EvEnd1 = WORMTRACK.end1(j,:);
    EvEnd2 = WORMTRACK.end2(j,:);
    
    centroids = [];
    
    end1gs = WORMTRACK.end1_g(j);
    end2gs = WORMTRACK.end2_g(j);
    
    ends_gs = [];
    
    if sum(isnan(EvEnd1))==0 && sum(isnan(EvEnd2))==0
        centroids = [EvEnd1; EvEnd2];
        ends_gs = [end1gs end2gs];
        %%% MAIN CODE FOR TRACKING HEAD AND TAIL
        predictNewLocationsOfTracks()
        [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment();
        updateAssignedTracks();
        updateUnassignedTracks();
        createNewTracks();
        %%%
        centroid_dist(j) = pdist([tracks(1).centroid(end,:);tracks(2).centroid(end,:)]);
        lastgood = 1;
        justbroke = 0;
    elseif ~isstruct(EvWorm) %if the worm was not segmented, put a 0 in the centroid_dist -- these are usually reorientations
        disp(['idx NOT SEGMENTED: ' num2str(j)]);
        centroid_dist(j) = 0;
    else
        error('not aware of this category');
    end
end

tracks_combined = [tracksThatLeft,tracks];
head = NaN(size(WORMTRACK.centroid));
tail = NaN(size(WORMTRACK.centroid));
head_gs = NaN(size(WORMTRACK.centroid,1),1);
tail_gs = NaN(size(WORMTRACK.centroid,1),1);
omega = zeros(size(WORMTRACK.centroid,1),1);

if isempty(tracks_combined) %if there were no segmented splines in this track, just return empty head and tail
    return
end

%filter out segmented frames in the midst of non-segmented frames -->
%centroid_dist is set to 0 at these points
nonzeroinds = find(centroid_dist>0);
x_diff = (diff(nonzeroinds)==1)';
f = find([false,x_diff]~=[x_diff,false]); %get consecutive series of frames nonzero
g = find(f(2:2:end)-f(1:2:end-1)<=3); %get indices of consecutive series less than or equal to 9 frames long
consec_starts = f(2*g-1);

all_frames_to_reset = [];
for r = 1:length(consec_starts)
    startframe = consec_starts(r);
    st = nonzeroinds(startframe);
    frames_to_reset = st;

    while centroid_dist(st+1)~=0
        frames_to_reset = [frames_to_reset; st+1];
        if st+1>=length(centroid_dist)
            break;
        else
            st = st+1;
        end
    end
    all_frames_to_reset = [all_frames_to_reset; frames_to_reset];
end

[cdh, cdh_ind] = hampel(centroid_dist,1);
cdh_ind = find(cdh_ind);
d = cdh(cdh_ind)==0; %these are the indices of centroid_dist that were replaced with zeros by the hampel filter.
all_frames_to_reset = unique([all_frames_to_reset; cdh_ind(d)]);
cd_copy = centroid_dist;
cd_copy(all_frames_to_reset)=0;
% cd_copy(isnan(cd_copy))=0;
centroid_dist = cd_copy;
% find split points
a = find(centroid_dist==0);
if ~isempty(a)
    split_starts = a((diff(a)~=1)');
    split_starts = unique([split_starts;a(end)]);%add back last zero member, sometimes lost
    a_diff = (diff(a)==1)';
    f = find([false,a_diff]~=[a_diff,false]);
    g = find(f(2:2:end)-f(1:2:end-1));
    split_ends = a(f(2*g-1));
    splits = unique(sort([1;split_starts;split_ends;length(centroid_dist)])); %these are all breakpoints
    split_int = [splits(1:end-1)+1 splits(2:end)-1];
    split_int(1) = 1; split_int(end) = length(centroid_dist); %make sure the first index is 1 and the last index is the last number in the array
    split_int(split_int(:,1)-split_int(:,2)>0,:) = []; %get rid of impossible intervals
    % and intervals to keep, those containing nonzero, nonnan entries
    keep_ints = [];
    for n = 1:size(split_int,1)
        c = centroid_dist(split_int(n,1):split_int(n,2));
        if sum(c==0)~=length(c)
            keep_ints = [keep_ints ; split_int(n,:)];
        end
    end
else %no need to split this track further
    keep_ints = [1 WORMTRACK.totalVisibleCount];
end

%%% FIND INTERVALS CONTAINING OMEGAS AMONG THE NON-SEGMENTED FRAMES
fr = tracks_combined(1).frame_rel;
seg_inds = [];
for o = 1:size(keep_ints,1)
    seg_inds = [seg_inds; (keep_ints(o,1):keep_ints(o,2))'];
end
non_seg_inds = setdiff((1:WORMTRACK.totalVisibleCount)',seg_inds);
c_diff = (diff(non_seg_inds)==1);
f = find([false;c_diff]~=[c_diff;false]);
g = find(f(2:2:end)-f(1:2:end-1));
non_seg_ints = [non_seg_inds(f(2*g-1)) non_seg_inds(f(2*g))]; %these are intervals, which may contain omegas

omega_frames = [];
if ~isempty(non_seg_ints)
    for l = 1:size(non_seg_ints,1)
        d = non_seg_ints(l,1):non_seg_ints(l,2);
        nonsegcount = 0;
        if isempty(d)
            continue;
        end
        ind_local = d;
        for m = 1:length(ind_local)
            wc = WORMTRACK.cropworm{ind_local(m)};
            wc_orig = WORMTRACK.cropworm_orig{ind_local(m)};
            [~,~,errNum,~] = segWorm_Elias(wc,wc_orig,ind_local(m), 0, 0);
            if errNum == 105
                nonsegcount = nonsegcount+1;
            end
        end
        if nonsegcount/length(d)>0.3
            omega_frames = [omega_frames; d'];
        end
    end
end

%%% STITCH HEAD AND TAIL SEGMENTS BACK TOGETHER %%%
fr = tracks_combined(1).frame_rel;
for o = 1:size(keep_ints,1)
    d = keep_ints(o,1):keep_ints(o,2);
    
    if isempty(d)
        continue;
    end
    %find the corresponding indices of d in tracks_combined.frame_rel
    [d_new,~,ind_local] = intersect(d,fr);
    d_new = d_new';
    
    %     if ind_local(end)>size(tracks_combined(1).centroid,1)
    %         ind_local = ind_local(1):size(tracks_combined(1).centroid,1);
    %     end
    
    cent1 = tracks_combined(1).centroid(ind_local,:);
    cent2 = tracks_combined(2).centroid(ind_local,:);
    
    end1_gs = tracks_combined(1).end_g(ind_local);
    end2_gs = tracks_combined(2).end_g(ind_local);
    
    a1 = squareform(pdist(cent1));
    b1 = a1(:,2:end);
    c1 = diag(b1);
    d1 = sum(c1);
    
    a2 = squareform(pdist(cent2));
    b2 = a2(:,2:end);
    c2 = diag(b2);
    d2 = sum(c2);
    
    if d2-d1>=5 % cent2 is the head
        head(d_new,:) = cent2;
        tail(d_new,:) = cent1;
        head_gs(d_new) = end2_gs;
        tail_gs(d_new) = end1_gs;
    elseif d1-d2>=5 % cent1 is the head
        head(d_new,:) = cent1;
        tail(d_new,:) = cent2;
        head_gs(d_new) = end1_gs;
        tail_gs(d_new) = end2_gs;
    else
        head(d_new,:) = NaN(size(d_new,2),2);
        tail(d_new,:) = NaN(size(d_new,2),2);
        head_gs(d_new) = NaN(size(d_new,2),1);
        tail_gs(d_new) = NaN(size(d_new,2),1);
        warning(['frame ' num2str(d_new(1)) 'to' num2str(d_new(end)) ': cannot discern H v T by displacement in this interval!']);
    end
end
%%%
%add a field called 'omega' where omega frames are marked - may update
%later with reversals, too

% omega(omega_frames) = cellstr(repmat('omega',size(omega_frames,1),1));
% a = find(contains(omega,'omega')) % this is how you find the indices with
% that behavioral annotation

omega(omega_frames) = 1;
omega = logical(omega);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [tracks, tracksThatLeft] =  initializeTracks()
        tracks = struct(...
            'id',{},...
            'x_off',{},...
            'y_off',{},...
            'centroid',{},...
            'end_g',{},...
            'kalmanFilter', {}, ...
            'age', {}, ...
            'totalVisibleCount', {}, ...
            'consecutiveInvisibleCount', {},...
            'frame_rel',{},...
            'frame_global',{});
        tracksThatLeft = struct(...
            'id',{},...
            'x_off',{},...
            'y_off',{},...
            'centroid',{},...
            'end_g',{},...
            'kalmanFilter', {}, ...
            'age', {}, ...
            'totalVisibleCount', {}, ...
            'consecutiveInvisibleCount', {},...
            'frame_rel',{},...
            'frame_global',{});
    end

    function predictNewLocationsOfTracks()
        for i = 1:length(tracks)
            % Predict the current location of the track.
            predict(tracks(i).kalmanFilter);
        end
    end

    function [assignments, unassignedTracks, unassignedDetections] = ...
            detectionToTrackAssignment()
        
        nTracks = length(tracks);
        nDetections = size(centroids, 1);
        
        % Compute the cost of assigning each detection to each track.
        cost = zeros(nTracks, nDetections);
        for i = 1:nTracks
            cost(i, :) = distance(tracks(i).kalmanFilter, centroids);
        end
        
        % Impose an additional constraint based on absolute euclidean distance --
        % set the cost function for any pair of detection and assignment to Inf if
        % the distance between the last centroid position and the detection is over
        % a certain distance threshold. the goal is to reduce some ambiguities.
%         lastpos = zeros(nTracks,2);
%         for i = 1:nTracks
%             lastpos(i,:) = tracks(i).centroid(end,:); %this is the kalman prediction
%         end
%         euc_dist = pdist2(lastpos,centroids,'euclidean');
%         [y, k1] = min(euc_dist);
%         k2 = (find(y>50));
%         if ~isempty(k2)
%             d_i = k2; d_j = k1(k2);
%             
%             ind2replace = sub2ind(size(cost),d_i,d_j);
%             cost(ind2replace) = Inf;
%         end
        
        
        % Solve the assignment problem.
        costOfNonAssignment = realmax/1000;
        if nDetections>2
            costOfNonAssignment = 50;
        end
        [assignments, unassignedTracks, unassignedDetections] = ...
            assignDetectionsToTracks(cost, costOfNonAssignment);
    end

    function updateAssignedTracks()
        numAssignedTracks = size(assignments, 1);
        for i = 1:numAssignedTracks
            trackIdx = assignments(i, 1);
            detectionIdx = assignments(i, 2);
            centroid = centroids(detectionIdx, :);
            end_g = ends_gs(detectionIdx);
            
            % Correct the estimate of the object's location
            % using the new detection.
            correct(tracks(trackIdx).kalmanFilter, centroid);
            
            centroid_upToNow = tracks(trackIdx).centroid;
            tracks(trackIdx).centroid = [centroid_upToNow; centroid];
            
            %Update the gs value associated with this endpoint
            gs_upToNow = tracks(trackIdx).end_g;
            tracks(trackIdx).end_g = [gs_upToNow; end_g];
            
            % Update track's age.
            tracks(trackIdx).age = tracks(trackIdx).age + 1;
            
            % Update visibility.
            tracks(trackIdx).totalVisibleCount = ...
                tracks(trackIdx).totalVisibleCount + 1;
            tracks(trackIdx).consecutiveInvisibleCount = 0;
            
            
            x_off_soFar = tracks(trackIdx).x_off;
            tracks(trackIdx).x_off = [x_off_soFar;x_off];
            
            y_off_soFar = tracks(trackIdx).y_off;
            tracks(trackIdx).y_off = [y_off_soFar;y_off];
            
            frame_rel_SoFar = tracks(trackIdx).frame_rel;
            tracks(trackIdx).frame_rel = [frame_rel_SoFar;j];
            
            frame_global_SoFar = tracks(trackIdx).frame_global;
            tracks(trackIdx).frame_global = [frame_global_SoFar; WORMTRACK.framesActive(j)];
        end
    end

    function updateUnassignedTracks()
        for i = 1:length(unassignedTracks)
            ind = unassignedTracks(i);
            
            tracks(ind).age = tracks(ind).age + 1;
            tracks(ind).consecutiveInvisibleCount = ...
                tracks(ind).consecutiveInvisibleCount + 1;

            x_off_soFar = tracks(ind).x_off;
            tracks(ind).x_off = [x_off_soFar;x_off];
            
            y_off_soFar = tracks(ind).y_off;
            tracks(ind).y_off = [y_off_soFar;y_off];
            
            frame_rel_SoFar = tracks(ind).frame_rel;
            tracks(ind).frame_rel = [frame_rel_SoFar;j];
            
            frame_global_SoFar = tracks(ind).frame_global;
            tracks(ind).frame_global = [frame_global_SoFar; WORMTRACK.framesActive(j)];
        end
        identifyTracksStillIn();
    end

    function identifyTracksStillIn()
        still_visible = NaN(length(tracks),1);
        for i = 1:length(tracks)
            still_visible(i) = ~(tracks(i).consecutiveInvisibleCount > 2); %was 2 before
        end
        stillGoodInds = find(still_visible);
        tracksStillGood = tracks(stillGoodInds); %these should be tracked
        newTracksThatLeft = tracks(setdiff(1:length(tracks),stillGoodInds));
        WorthSavingTracksThatLeft_idx = [newTracksThatLeft(:).age]'>0; %only save tracks that left which were at least X frames long
        newTracksThatLeft = newTracksThatLeft(WorthSavingTracksThatLeft_idx);
        tracksThatLeft = [tracksThatLeft, newTracksThatLeft];
        tracks = tracksStillGood;
    end

    function createNewTracks()
        if length(unassignedDetections)>2
            unassignedDetections = unassignedDetections(1:2);
        end
        if isempty(tracks)
            centroids_loc = centroids(unassignedDetections, :);
            ends_g_loc = ends_gs(unassignedDetections);
        else
            centroids_loc = [];
            ends_g_loc = [];
        end
        
        for i = 1:size(centroids_loc, 1)
            
            centroid = centroids_loc(i,:);
            end_g = ends_g_loc(i);
            
            % Create a Kalman filter object.
            kalmanFilter = configureKalmanFilter('ConstantVelocity', ...
                centroid, [200, 50], [100, 25], 100);
            
            % Create a new track.
            newTrack = struct(...
                'id', nextId, ...
                'x_off',x_off,...
                'y_off',y_off,...
                'centroid', centroid, ...
                'end_g',end_g,...
                'kalmanFilter', kalmanFilter, ...
                'age', 1, ...
                'totalVisibleCount', 1, ...
                'consecutiveInvisibleCount', 0, ...
                'frame_rel',j,...
                'frame_global',WORMTRACK.framesActive(j));
            
            % Add it to the array of tracks.
            tracks(end + 1) = newTrack;
            
            % Increment the next id.
            nextId = nextId + 1;
        end
    end

    function displayTrackingResults()
        minVisibleCount = 0;
        orig = getWormCrop();
        if ~isempty(tracks)
            
            % Display the objects. If an object has not been detected
            % in this frame, display its predicted bounding box.
            
            % Get positions.
            for k = 1:length(tracks)
                tmp_centroid = tracks(k).centroid;
                curr_centroids = tmp_centroid(size(tmp_centroid,1),:);
                curr_centroids = [curr_centroids(1)-x_off curr_centroids(2)-y_off];%put back in worm coordinates for display
                centroids_to_draw = [curr_centroids 3*ones(size(curr_centroids,1),1)];
                
                % Get ids.
                ids = int32([tracks(k).id]);
                
                % Create labels for objects indicating the ones for
                % which we display the predicted rather than the actual
                % location.
                labels = cellstr(int2str(ids'));
                predictedTrackInds = ...
                    [tracks(k).consecutiveInvisibleCount] > 0;
                isPredicted = cell(size(labels));
                isPredicted(predictedTrackInds) = {' predicted'};
                labels = strcat(labels, isPredicted);
                
                % Draw the objects on the frame.
                orig = insertShape(orig,'FilledCircle',centroids_to_draw,'LineWidth',1,'Color','magenta');
                orig = insertObjectAnnotation(orig, 'Circle', ...
                    centroids_to_draw, labels,'FontSize',8);
            end
            
        end
        
        % Display the mask and the frame.
        origPlayer.step(orig);
        release(origPlayer);
        pause();
        %         showtracking = 0;
    end

    function orig = getWormCrop()
        tmpspline = WORMTRACK.spline{j}';
        tmpspline(1,:) = tmpspline(1,:)-x_off;
        tmpspline(2,:) = tmpspline(2,:)-y_off;
        %express splines as an M-by-2L matrix, where each row is a vector representing a polyline with L number of vertices.
        t(1:2:length(tmpspline)*2) = tmpspline(1,:);
        t(2:2:length(tmpspline)*2) = tmpspline(2,:);
        
        evspline = EvSkel';
        x(1:2:length(evspline)*2) = evspline(1,:);
        x(2:2:length(evspline)*2) = evspline(2,:);
        
        orig = WORMTRACK.cropworm_orig{j};
        orig = insertShape(orig,'Line',t,'LineWidth',1,'Color','cyan');
        if sum(sum(isnan(x)))==0
            orig = insertShape(orig,'Line',x,'LineWidth',1,'Color','magenta');
        end
    end

%     function splitTracks()
%         if lastgood %if the last frame had a spline and this one doesn't
%             badcount = 1;
%         else
%             badcount = badcount+1;
%             if badcount>2 && ~justbroke %after 2 NaNs or no spline, split tracks, just to give buffer for a blip
%                 tracksThatLeft = [tracksThatLeft, tracks];
%                 tracks = tracks([]);
%                 justbroke = 1; %this get's reset as soon as there is a good observation
%             end
%         end
%         lastgood = 0;
%     end

end










