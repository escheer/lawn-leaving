%% Display Tracking Results
% The |displayTrackingResults| function draws a bounding box and label ID
% for each track on the video frame and the foreground mask. It then
% displays the frame and the mask in their respective video players.

% Elias Scheer
% 1-12-17
% borrowed from MultiObjectTracking example

function displayTrackingResults_wholevideo_lawn(blobPlayer, bboxes, trackPlayer, tracks, thresh_cleaned, ev_ho_x, ev_ho_y)%, image_y_dim, curr_enter_events, curr_exit_events)
%blob player
blobs = insertShape(double(thresh_cleaned), 'Rectangle', bboxes, 'Color', 'green','LineWidth',2);
blobPlayer.step(blobs);

masked_worms = uint8(repmat(thresh_cleaned, [1, 1, 3])) .* 255;

if ~isempty(tracks)
    
    % Get centroids and bounding boxes.
    curr_centroids = NaN(length(tracks),2);
    bboxes = NaN(length(tracks),4);
    for j = 1:length(tracks)
        tmp_centroid = tracks(j).centroid;
        curr_centroids(j,:) = tmp_centroid(size(tmp_centroid,1),:);
        curr_bbox = tracks(j).bbox;
        bboxes(j,:) = curr_bbox(size(curr_bbox,1),:);
    end
    % Get ids.
    ids = int32([tracks(:).id]);
    
    % Create labels for objects indicating the ones for
    % which we display the predicted rather than the actual
    % location.
    labels = cellstr(int2str(ids'));
    predictedTrackInds = ...
        [tracks(:).consecutiveInvisibleCount] > 0;
    isPredicted = cell(size(labels));
    isPredicted(predictedTrackInds) = {' predicted'};
    labels = strcat(labels, isPredicted);
    
    % write ENT or EXT on entering or exiting tracks, respectively
%     entering = cell(size(labels));
%     idx = ismember(ids,curr_enter_events); %compare the track ids to find the correct indices
%     c = 1:length(ids);
%     d_ent = c(idx);
%     entering(d_ent) = {' ENTER'};
%     labels = strcat(labels, entering);
    
%     exiting = cell(size(labels));
%     idx = ismember(ids,curr_exit_events);
%     c = 1:length(ids);
%     d_ext = c(idx);
%     exiting(d_ext) = {' EXIT'};
%     labels = strcat(labels, exiting);
%     
    % Draw the objects on the mask.
    centroids_to_draw = [curr_centroids 3*ones(size(curr_centroids,1),1)];
    masked_worms = insertShape(masked_worms,'FilledCircle',centroids_to_draw,'LineWidth',1,'Color','magenta');
    for j = 1:length(tracks) %plots the entire trailing centroid path before now
        masked_worms = insertMarker(masked_worms,tracks(j).centroid,'*','color','blue','size',1);
    end
    masked_worms = insertObjectAnnotation(masked_worms, 'rectangle', bboxes, labels,'Color','magenta');
%     if ~isempty(curr_enter_events)
%         masked_worms = insertObjectAnnotation(masked_worms, 'rectangle', bboxes(d_ent,:), labels(d_ent),'Color','green');
%     end
%     if ~isempty(curr_exit_events)
%         masked_worms = insertObjectAnnotation(masked_worms, 'rectangle', bboxes(d_ext,:), labels(d_ext),'Color','red');
%     end
    masked_worms = insertText(masked_worms, [10 10], length(tracks), 'BoxOpacity', 1, 'FontSize', 45);
%     masked_worms = insertMarker(masked_worms,[image_y_dim - ev_ho_y' ev_ho_x'],'+','color','red','size',1);
    masked_worms = insertMarker(masked_worms,[ev_ho_x ev_ho_y],'+','color','red','size',1);
end

% Display the mask and the frame.
trackPlayer.step(masked_worms);

end