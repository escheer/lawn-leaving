%% Display Tracking Results
% The |displayTrackingResults| function draws a bounding box and label ID
% for each track on the video frame and the foreground mask. It then
% displays the frame and the mask in their respective video players.

% Elias Scheer
% 1-12-17
% borrowed from MultiObjectTracking example

function displayTrackingResults_postproc(tracks,trackrange,framerange)
close all;
if nargin<3
    doall = 1;
else
    doall = 0;
end

origPlayer = vision.VideoPlayer('Position', [120, 160, 100, 100]);
threshPlayer = vision.VideoPlayer('Position', [970, 160, 100, 100]);

if doall
    trackrange = 1:length(tracks);
end
for i = trackrange
%     disp(['i: ' num2str(i)]);
    curr_track = tracks(i);
    if doall
        framerange = 1:curr_track.age;
    end
    for j = framerange
        disp(['j: ' num2str(j)]);
        orig = curr_track.cropworm_orig{j};
        thresh = im2uint8(curr_track.cropworm{j});
        
        frame = curr_track.framesActive(j);
        [worm errNum errMsg] = segWorm(orig, frame, 0, 1);
        
        x_off = double(curr_track.bbox(j,1));
        y_off = double(curr_track.bbox(j,2));
        head = curr_track.head;
        tail = curr_track.tail;
        
        tmpspline = curr_track.spline{j}';
        if ~isempty(tmpspline) && sum(sum(isnan(tmpspline)))<1
            tmpspline(1,:) = tmpspline(1,:)-x_off;
            tmpspline(2,:) = tmpspline(2,:)-y_off;
            %express splines as an M-by-2L matrix, where each row is a vector representing a polyline with L number of vertices.
            t(1:2:length(tmpspline)*2) = tmpspline(1,:);
            t(2:2:length(tmpspline)*2) = tmpspline(2,:);

            if ~isnan(head) && ~isnan(tail)
                head_x = double(curr_track.head(j,1)-x_off);
                head_y = double(curr_track.head(j,2)-y_off);
                tail_x = double(curr_track.tail(j,1)-x_off);
                tail_y = double(curr_track.tail(j,2)-y_off);
                
                orig = insertMarker(orig,[head_x head_y],'+','color','green','size',1);
                orig = insertMarker(orig,[tail_x tail_y],'+','color','red','size',1);
                thresh = insertMarker(thresh,[head_x head_y],'+','color','green','size',2);
                thresh = insertMarker(thresh,[tail_x tail_y],'+','color','red','size',2);
            end
            
                orig = insertShape(orig,'Line',t,'LineWidth',1,'Color','cyan');
                origPlayer.step(orig);
                release(origPlayer);
                
                thresh = insertShape(thresh,'Line',t,'LineWidth',1,'Color','cyan');
                threshPlayer.step(thresh);
                release(threshPlayer);
      
        end
        pause();
        close all;
    end
end


end
