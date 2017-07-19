%% Update Assigned Tracks
% The |updateAssignedTracks| function updates each assigned track with the
% corresponding detection. It calls the |correct| method of
% |vision.KalmanFilter| to correct the location estimate. Next, it stores
% the new bounding box, and increases the age of the track and the total
% visible count by 1. Finally, the function sets the invisible count to 0.

% Elias Scheer
% 1-12-17
% borrowed from MultiObjectTracking example

function tracks = updateAssignedTracks_wholevideo(tracks,assignments,centroids,bboxes,curr_frame)
numAssignedTracks = size(assignments, 1);
for i = 1:numAssignedTracks
    trackIdx = assignments(i, 1);
    detectionIdx = assignments(i, 2);
    centroid = centroids(detectionIdx, :);
    bbox = bboxes(detectionIdx, :);
    
    % Correct the estimate of the object's location
    % using the new detection.
    correct(tracks(trackIdx).kalmanFilter, centroid);
    
    % Replace predicted bounding box with detected
    % bounding box & append the new detection to the past observations
    bbox_upToNow = tracks(trackIdx).bbox;
    bbox_upToNow(size(bbox_upToNow,1),:) = bbox; %remove predicted one and correct it with new one
    tracks(trackIdx).bbox = bbox_upToNow;
    centroid_upToNow = tracks(trackIdx).centroid;
    centroid_upToNow(size(bbox_upToNow,1),:) = centroid; %remove predicted one and correct it with new one
    tracks(trackIdx).centroid = centroid_upToNow;
    
    % Update track's age.
    tracks(trackIdx).age = tracks(trackIdx).age + 1;
    
    % Update visibility.
    tracks(trackIdx).totalVisibleCount = ...
        tracks(trackIdx).totalVisibleCount + 1;
    tracks(trackIdx).consecutiveInvisibleCount = 0;
    
    %Update frames active.
    tracks(trackIdx).framesActive = [tracks(trackIdx).framesActive curr_frame]; 
    
end
end