function [enter_events, exit_events] = countCrossingTracks_after( tracks,tracksThatLeft, ev_ho_x, ev_ho_y)
%countCrossingTracks_after.m This function takes in tracks and
%tracksThatLeft and looks for tracks which intersect with the event
%horizon. Entering or exiting events are marked as the frame when the entire bounding box changes side.
%It then classifies these tracks as entering or exiting based on
%the direction of the centroid path.
% Elias Scheer 04-24-17
%% 0. get data in useful format
x = ev_ho_x;
x_segs = [[x(1:end-1) x(2:end)];[x(end) x(1)]]; %turn collection of vertices into line segments
y = ev_ho_y;
y_segs = [[y(1:end-1) y(2:end)];[y(end) y(1)]];
linesegs = [x_segs(:,1) y_segs(:,1) x_segs(:,2) y_segs(:,2)];

tracksToCheck = [tracks, tracksThatLeft]; %merge tracks and tracksThatLeft
%% 1.Identify tracks crossing the event horizon.
% enter_fig = figure();hold on;plot(x,y);axis equal;title('ENTERING');
% exit_fig = figure();hold on;plot(x,y);axis equal;title('EXITING');

centroids = {tracksToCheck(:).centroid};
enter_events = []; %this will be filled in by [framenum trackid]
exit_events = []; %same

for i = 1:length(tracksToCheck)
    %     if i == find([tracksToCheck.id] == 3217)
    %         disp('debug');
    %     end
    trackid = tracksToCheck(i).id;
    firstframe = tracksToCheck(i).framesActive(1);
    lastframe = tracksToCheck(i).framesActive(end);
   
    tmp = double(centroids{i}); %convert to double for all subsequent tests
    centroid_segs = [tmp(1:end-1,:) tmp(2:end,:)];
    out = lineSegmentIntersect(linesegs,centroid_segs);
    num_int = sum(sum(out.intAdjacencyMatrix));
    if num_int == 0 %if this track had no crossings, go to the next
        continue;
    else
        adjmat = out.intAdjacencyMatrix;
        [~,seg_idx] = find(~isnan(adjmat)&adjmat-0~=0); %this is the index of the segment of the centroid path which crossed
        crossing_frames = seg_idx+firstframe-1; %need to subtract 1 to align everything
        
        % now check for which of these crossing events result in full
        % crossing of bbox
        
        %check that the track exists for at least X frames before and after
        %the first and last crossings respectively
        %         framesToCheck = min(crossing_frames)-5:max(crossing_frames)+5;
        
        framesToCheck = firstframe:lastframe;
        idxToCheck = framesToCheck-firstframe+1; %track indices
        if ismember(framesToCheck(1),tracksToCheck(i).framesActive)&&ismember(framesToCheck(end),tracksToCheck(i).framesActive)
            bbox = tracksToCheck(i).bbox(idxToCheck,:);
            bbox_verts =   [[bbox(:,1) bbox(:,2)]...          % top left
                [bbox(:,1) bbox(:,2)+bbox(:,4)]...            % bottom left
                [bbox(:,1)+bbox(:,3) bbox(:,2)]...            % top right
                [bbox(:,1)+bbox(:,3) bbox(:,2)+bbox(:,4)]];   % bottom right
            
            topleft_out = ~inpolygon(double(bbox_verts(:,1)),double(bbox_verts(:,2)),x,y);
            bottomleft_out = ~inpolygon(double(bbox_verts(:,3)),double(bbox_verts(:,4)),x,y);
            topright_out = ~inpolygon(double(bbox_verts(:,5)),double(bbox_verts(:,6)),x,y);
            bottomright_out = ~inpolygon(double(bbox_verts(:,7)),double(bbox_verts(:,8)),x,y);
            
            bboxes_out = topleft_out&bottomleft_out&topright_out&bottomright_out; %subtrack indices
            bboxes_in = ~topleft_out&~bottomleft_out&~topright_out&~bottomright_out; %subtrack indices
            
            inside_to_X = find([0 diff(bboxes_in)']==-1)'+idxToCheck(1)-1; %put it back in track indices
            X_to_inside = find([0 diff(bboxes_in)']==1)'+idxToCheck(1)-1;
            outside_to_X = find([0 diff(bboxes_out)']==-1)'+idxToCheck(1)-1;
            X_to_outside = find([0 diff(bboxes_out)']==1)'+idxToCheck(1)-1;
            
            % by convention, we refer to leaving events as +1 and entering as -1
            glom =  [[inside_to_X ones(size(inside_to_X))];...
                [X_to_outside ones(size(X_to_outside))];
                [outside_to_X -1*ones(size(outside_to_X))];
                [X_to_inside -1*ones(size(X_to_inside))]];
            glom = sortrows(glom,1);
            consec = [1;diff(glom(:,2))]==0;
            cross_inds = [glom(consec,1) glom(consec,2)]; %look for consecutive upsteps and downsteps
            
            enter_ind = cross_inds(cross_inds(:,2)==-1,1);%two consecutive downsteps = entering
            exit_ind = cross_inds(cross_inds(:,2)==1,1); %two consecutive upsteps = exiting
            
            if isempty(enter_ind)&&isempty(exit_ind) %if none of these centroid crossings amounted to a bbox crossing, skip track
                continue;
            end
            if ~isempty(enter_ind) %we had some entering events, add them to the list
                enter_frame = enter_ind+firstframe-1; %put it back in frame indices
                newEnteringEvents = [enter_frame trackid.*ones(size(enter_frame))];
                enter_events = [enter_events ; newEnteringEvents];
                %                 figure(enter_fig);
                %                 plot(tracksToCheck(i).centroid(:,1),tracksToCheck(i).centroid(:,2));
                %                 scatter(tracksToCheck(i).centroid(1,1),tracksToCheck(i).centroid(1,2),'g');
                %                 scatter(tracksToCheck(i).centroid(end,1),tracksToCheck(i).centroid(end,2),'r');
                %                 enter_bbox = tracksToCheck(i).bbox(enter_ind,:);
                %                 for k = 1:size(enter_bbox,1)
                %                     rectangle('Position',enter_bbox(k,:));
                %                     text(double(enter_bbox(k,1)),double(enter_bbox(k,2)),[num2str(enter_frame) ', ' num2str(trackid)]);
                %                     pause();
                %                 end
            end
            if ~isempty(exit_ind) %we had some exiting events, add them to the list
                exit_frame = exit_ind+firstframe-1;
                newExitingEvents = [exit_frame trackid.*ones(size(exit_frame))];
                exit_events = [exit_events ; newExitingEvents];
                %                 figure(exit_fig);
                %                 plot(tracksToCheck(i).centroid(:,1),tracksToCheck(i).centroid(:,2));
                %                 scatter(tracksToCheck(i).centroid(1,1),tracksToCheck(i).centroid(1,2),'g');
                %                 scatter(tracksToCheck(i).centroid(end,1),tracksToCheck(i).centroid(end,2),'r');
                %                 exit_bbox = tracksToCheck(i).bbox(exit_ind,:);
                %                 for k = 1:size(exit_bbox,1)
                %                     rectangle('Position',exit_bbox(k,:));
                %                     text(double(exit_bbox(k,1)),double(exit_bbox(k,2)),[num2str(exit_frame) ', ' num2str(trackid)]);
                %                     pause();
                %                 end
            end
        else
            continue; %track doesn't have enough buffer on either side of crossing events to evaluate whether there was a crossing
        end
    end
end

close all;
end

