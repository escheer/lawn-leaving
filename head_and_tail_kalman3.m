function [omega, head, tail, head_gs, tail_gs] = head_and_tail_kalman3( WORMTRACK )
%head_and_tail_kalman.m This function tracks each endpoint of the spline
%independently with kalman filters.

origPlayer = vision.VideoPlayer('Position', [120, 160, 100, 100]);

[tracks,tracksThatLeft] =  initializeTracks();
nextId = 1;

nonseg = zeros(WORMTRACK.totalVisibleCount,1);
for j = 1:WORMTRACK.totalVisibleCount
    if j == 3343
        disp('debug');
    end
    disp(['idx: ' num2str(j)]);

    x_off = double(WORMTRACK.bbox(j,1))-5; %x and y offsets of the bounding box
    y_off = double(WORMTRACK.bbox(j,2))-5;
    
    EvWorm = WORMTRACK.worm{j};
    
    EvEnd1 = WORMTRACK.end1(j,:); %one end of the thresholded worm, and the other
    EvEnd2 = WORMTRACK.end2(j,:);
    
    centroids = [];
    
    end1gs = WORMTRACK.end1_g(j); %grayscale of either end of the worm
    end2gs = WORMTRACK.end2_g(j);
    
    ends_gs = [];
    
    if sum(isnan(EvEnd1))==0 && sum(isnan(EvEnd2))==0
        centroids = [EvEnd1; EvEnd2];
        ends_gs = [end1gs end2gs];
        %%% MAIN CODE FOR TRACKING HEAD AND TAIL
        predictNewLocationsOfTracks();
        [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment();
        updateAssignedTracks();
        updateUnassignedTracks();
        updateTrackPartner();% update track partners before the creation of new tracks
        createNewTracks();
        
    elseif ~isstruct(EvWorm) %if the worm was not segmented, put a 0 in the centroid_dist -- these are usually reorientations
        disp(['idx NOT SEGMENTED: ' num2str(j)]);
        nonseg(j) = 1;
        %end tracks here
        if ~isempty(tracks)
            tracksThatLeft = [tracksThatLeft, tracks]; %#ok<AGROW>
            [tracks,~] =  initializeTracks();
        end
    else
        error('not aware of this category: either worm has endpoints or the EvWorm struct should be empty');
    end
end

tracks_combined = [tracksThatLeft,tracks]; %combine the tracked trajectories of endpoints
tracks_combined = tracks_combined([tracks_combined(:).age]'>=9); %only tracks at least 9 frames long

omega = zeros(size(WORMTRACK.centroid,1),1); %pre-allocate
head = NaN(size(WORMTRACK.centroid));      
tail = NaN(size(WORMTRACK.centroid));
head_gs = NaN(size(WORMTRACK.centroid,1),1);
tail_gs = NaN(size(WORMTRACK.centroid,1),1);

if isempty(tracks_combined) %if there were no segmented splines in this track, just return empty head and tail
    return
end

%%% FIND INTERVALS CONTAINING OMEGAS AMONG THE NON-SEGMENTED FRAMES
seg_inds = find(nonseg==0);
%fill in any missing frames to complete the intervals
x_diff = (diff(seg_inds)==1);
f = find([false;x_diff]~=[x_diff;false]);
g = find(f(2:2:end)-f(1:2:end-1)<=3);
consec_starts = f(2*g-1);
all_frames_to_reset = [];
for r = 1:length(consec_starts)
    startframe = consec_starts(r);
    st = seg_inds(startframe);
    frames_to_reset = st;

    while nonseg(st+1)~=1
        frames_to_reset = [frames_to_reset; st+1];
        if st+1>=length(nonseg)
            break;
        else
            st = st+1;
        end
    end
    all_frames_to_reset = [all_frames_to_reset; frames_to_reset];
end
[ns, ns_ind] = hampel(nonseg,2);
ns_ind = find(ns_ind);
d = ns(ns_ind)==1; %these are the indices of centroid_dist that were replaced with zeros by the hampel filter.
all_frames_to_reset = unique([all_frames_to_reset; ns_ind(d)]);
nonseg(all_frames_to_reset)=1;

%now create intervals
non_seg_inds = find(nonseg);
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
omega(omega_frames) = 1;
omega = logical(omega);

%%% STITCH HEAD AND TAIL SEGMENTS BACK TOGETHER %%%
id = [tracks_combined(:).id]';
id = [id zeros(size(id))]; %put a truth table next to 
for trk = 1:length(tracks_combined)
%     disp(['track ID is ' num2str(trk)]);
    curr_track = tracks_combined(trk);
    pID = curr_track.partnerID;
    allsame = sum(pID==pID(1))==length(pID); %check that the track only has 1 partner. this should always be true
    if allsame
        pID = pID(1);
        if id(id(:,1)==pID,2)==0 %if this track has not yet been processed
%             disp(['other track ID is ' num2str(pID)]);
            other_track = tracks_combined(id==pID);
            fr_curr = curr_track.frame_rel;
            fr_other = curr_track.frame_rel;
            if isequal(fr_curr, fr_other)
                fr = fr_curr;
            else
                error('track and partner track should exist for the same frames!');
            end
            cent1 = curr_track.centroid;
            cent2 = other_track.centroid;
            
            end1_gs = curr_track.end_g;
            end2_gs = other_track.end_g;
            
            a1 = squareform(pdist(cent1));
            b1 = a1(:,2:end);
            c1 = diag(b1);
            d1 = nansum(c1);
            
            a2 = squareform(pdist(cent2));
            b2 = a2(:,2:end);
            c2 = diag(b2);
            d2 = nansum(c2);
            
            if d2-d1>=5 % cent2 is the head
                head(fr,:) = cent2;
                tail(fr,:) = cent1;
                head_gs(fr) = end2_gs;
                tail_gs(fr) = end1_gs;
            elseif d1-d2>=5 % cent1 is the head
                head(fr,:) = cent1;
                tail(fr,:) = cent2;
                head_gs(fr) = end1_gs;
                tail_gs(fr) = end2_gs;
            else
                head(fr,:) = NaN(size(fr,1),2);
                tail(fr,:) = NaN(size(fr,1),2);
                head_gs(fr) = NaN(size(fr,1),1);
                tail_gs(fr) = NaN(size(fr,1),1);
                warning(['frame ' num2str(fr(1)) 'to' num2str(fr(end)) ': cannot discern H v T by displacement in this interval!']);
            end
            
            %mark each track as having been processed
            id(id==curr_track.id,2)=1;
            id(id==pID,2)=1;
        else
%             disp('skip this one, already processed...');
            continue;
        end
    else
        error('each track should always have only one partner!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [tracks, tracksThatLeft] =  initializeTracks()
        tracks = struct(...
            'id',{},...
            'partnerID',{},...
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
            'partnerID',{},...
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

    function updateTrackPartner()
        nTracks = length(tracks);
        if nTracks == 2 %if there are two tracks, partner them
            tracks(1).partnerID = [tracks(1).partnerID ; tracks(2).id];
            tracks(2).partnerID = [tracks(2).partnerID ; tracks(1).id];
        else
            for i = 1:nTracks %otherwise, the track has no partner on this frame
                tracks(i).partnerID = [tracks(i).partnerID ; NaN];
            end
        end
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
        
        % Remove any centroid detections that result in minimum
        % displacement of X (thresh for illegal tracking addition)
        lastpos = zeros(nTracks,2);
        for i = 1:nTracks
            lastpos(i,:) = tracks(i).centroid(end,:); %this is the kalman prediction
        end
        euc_dist = pdist2(lastpos,centroids,'euclidean');
        [y, i] = min(euc_dist);
        deletethis = false;
        k2 = (find(y>20)); %distance threshold 20
        if ~isempty(k2)
            deletethis = true;
        end
        if ~isempty(i)
            if i(1)==i(2)
                deletethis = true;
            end
        end
        
        if nTracks==0 && size(centroids, 1)==2 % if this is the first detection of a new track, don't allow detections to be too close together -- these are mistakes
            htdist = pdist2(centroids(1,:),centroids(2,:));
            if htdist < 35 %hard-coded threshold is far below worm length in both shays and my setups
                deletethis = true;
            end
        end
        
        if deletethis %remove centroid detections that result in a minimum displacement of X, also remove those detections in which the same last point is closer to both new detections
%             centroids(k2,:) = [];
            centroids = [];
        end
        
        nDetections = size(centroids, 1);
        
        % Compute the cost of assigning each detection to each track.
        cost = zeros(nTracks, nDetections);
        for i = 1:nTracks
            cost(i, :) = distance(tracks(i).kalmanFilter, centroids);
        end
        
        % Solve the assignment problem.
        costOfNonAssignment = realmax/1000;
        if nDetections~=2
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
            
            % Update the centroid positions
            centroid_upToNow = tracks(trackIdx).centroid;
            tracks(trackIdx).centroid = [centroid_upToNow; centroid];
            
            % Update the gs value associated with this endpoint
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
            
%             tracks(ind).age = tracks(ind).age + 1;
            tracks(ind).consecutiveInvisibleCount = ...
                tracks(ind).consecutiveInvisibleCount + 1;
            
%             % Update the centroid position with Nan
%             centroid_upToNow = tracks(ind).centroid;
%             tracks(ind).centroid = [centroid_upToNow; NaN NaN];
%             
%             %Update the gs value associated with this endpoint with NaN
%             gs_upToNow = tracks(ind).end_g;
%             tracks(ind).end_g = [gs_upToNow; NaN];
% 
%             x_off_soFar = tracks(ind).x_off;
%             tracks(ind).x_off = [x_off_soFar;x_off];
%             
%             y_off_soFar = tracks(ind).y_off;
%             tracks(ind).y_off = [y_off_soFar;y_off];
%             
%             frame_rel_SoFar = tracks(ind).frame_rel;
%             tracks(ind).frame_rel = [frame_rel_SoFar;j];
%             
%             frame_global_SoFar = tracks(ind).frame_global;
%             tracks(ind).frame_global = [frame_global_SoFar; WORMTRACK.framesActive(j)];
        end
        identifyTracksStillIn();
    end

    function identifyTracksStillIn()
        still_visible = NaN(length(tracks),1);
        for i = 1:length(tracks)
            still_visible(i) = ~(tracks(i).consecutiveInvisibleCount > 0); %was 2 before
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
        if length(tracks)<2 %is this right?
            centroids_loc = centroids(unassignedDetections, :);
            ends_g_loc = ends_gs(unassignedDetections);
        else
            centroids_loc = [];
            ends_g_loc = [];
        end
        
        partnertrack_id = [];
        if size(centroids_loc,1)==1 %we are making 1 new track, keep track of the existing one
            partnertrack_id = tracks.id;
        elseif size(centroids_loc,1)==2 %we are making 2 new tracks
            partnertrack_id = [nextId+1; nextId];
        else
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
                'partnerID',partnertrack_id(i),...
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

end










