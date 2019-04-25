function playframes( track, bg_struct, indices, coherence_thresh, speed_thresh, behavFLAG, deployable )
%PLAYFRAMES.M This little function is a wrapper for videoplayer to watch
%defined set of frames using only the cropped images of the worm after
%tracking is over.
% frames is in track indices -- this will be translated back to video frame
% indices in this function

if ischar(indices)
    if strcmp(indices,'all')
        indices = 1:track.age;
    else
        error('if string is specified, must be "all" ');
    end
end

if ~deployable
    whole_player = vision.VideoPlayer('Position', [1043, 32, 816, 593]);
else
    whole_player = vision.DeployableVideoPlayer('Size','Custom','CustomSize',[700 700]);
end

[head_vector, cent_vector, tail_vector, speed_smooth, moving_forward, moving_backward] = getforwardreverse( track.centroid,track.head,track.tail,track.speed, coherence_thresh, speed_thresh );

if behavFLAG
    head_pokes = track.head_pokes;
    if isempty(head_pokes)
        head_pokes = zeros(length(track.framesActive),1);
    end
    head_pokes_rev = track.head_poke_reversal; %changed 04/30/2018
    if isempty(head_pokes_rev)
        head_pokes_rev = zeros(length(track.framesActive),1);
    end
    lawn_entries = track.lawn_entries;
    lawn_exits = track.lawn_exits;
    roam_dwell = track.roamdwell_hmm;
end
currvid = track.videoname{indices(1)};
v = VideoReader(currvid);

last_x_offset = bg_struct(track.bgvidindex(indices(1))).region_rounded(1);
last_y_offset = bg_struct(track.bgvidindex(indices(1))).region_rounded(2);

for i = 1:length(indices)
    os_wait = 0;
    wait = 0;
    curr_idx = indices(i);
    
    disp(curr_idx);
    
    videoname = track.videoname{curr_idx};
    videoframenum = track.videoframe(curr_idx);
    if ~strcmp(currvid,videoname)
        clear v;
        v = VideoReader(videoname); %#ok<TNMLP>
    end
    wholeframe = read(v,videoframenum);
    wholeframe = insertText(wholeframe, [10 10], curr_idx, 'BoxOpacity', 1, 'FontSize', 45);
    
    %02/21/2018 there appears to be an off by one bug in bgvidindex -- not
    %obvious where it comes from but just look to the next one for now as
    %long as not at the end.
    
    %     if curr_idx<length(track.bgvidindex)
    %         next = curr_idx+1;
    %     else
    %         next = curr_idx;
    %     end
    next = curr_idx;
    
    x_offset = bg_struct(track.bgvidindex(next)).region_rounded(1);
    y_offset = bg_struct(track.bgvidindex(next)).region_rounded(2);
    
    if x_offset < last_x_offset
        diff = abs(x_offset - last_x_offset);
        wholeframe = insertText(wholeframe,[400 100],['shift LEFT ' num2str(diff)],'FontSize',45,'BoxColor','white','BoxOpacity',1);
        os_wait = 1;
    elseif x_offset > last_x_offset
        diff = abs(x_offset - last_x_offset);
        wholeframe = insertText(wholeframe,[400 100],['shift RIGHT ' num2str(diff)],'FontSize',45,'BoxColor','white','BoxOpacity',1);
        os_wait = 1;
    end
    if y_offset < last_y_offset
        diff = abs(y_offset - last_y_offset);
        wholeframe = insertText(wholeframe,[400 200],['shift UP ' num2str(diff)],'FontSize',45,'BoxColor','white','BoxOpacity',1);
        os_wait = 1;
    elseif y_offset > last_y_offset
        diff = abs(y_offset - last_y_offset);
        wholeframe = insertText(wholeframe,[400 200],['shift DOWN ' num2str(diff)],'FontSize',45,'BoxColor','white','BoxOpacity',1);
        os_wait = 1;
    end
    
    last_x_offset = x_offset;
    last_y_offset = y_offset;
    
    % plot event horizon
    eh = bg_struct(track.bgvidindex(next)).ev_ho;
    ev_ho = zeros(1,2*length(eh));
    ev_ho(1:2:length(ev_ho))=eh(:,1);
    ev_ho(2:2:length(ev_ho))=eh(:,2);
    wholeframe = insertShape(wholeframe,'Line',ev_ho,'LineWidth',2,'Color','w');
    
    %plot centroid and speed
    centroid = track.centroid(curr_idx,:);
    centroid = [centroid(:,1)+x_offset centroid(:,2)+y_offset];
    wholeframe = insertShape(wholeframe,'FilledCircle',[centroid 5],'LineWidth',2,'Color','magenta');
    
    % insert head, centroid, and tail movement vectors
    vizscale = 5;
    curr_cv = cent_vector(curr_idx,:).*vizscale;
    if sum(isnan(curr_cv))==0
        %         wholeframe = insertShape(wholeframe,'Line',[centroid centroid+curr_cv],'LineWidth',2,'Color','magenta');
    end
    
    %UNDO THIS! 01/14/2019
% %     head = track.head(curr_idx,:);
% %     if ~isempty(head) && sum(isnan(head))==0
% %         head = [head(:,1)+x_offset head(:,2)+y_offset]; 
% %         wholeframe = insertShape(wholeframe,'FilledCircle',[head 3],'LineWidth',2,'Color','green');
% %         curr_hv = head_vector(curr_idx,:).*vizscale;
% %         if sum(isnan(curr_hv))==0
% %             %             wholeframe = insertShape(wholeframe,'Line',[head head+curr_hv],'LineWidth',2,'Color','green');
% %         end
% %     end
% %     tail = track.tail(curr_idx,:);
% %     if ~isempty(tail) && sum(isnan(tail))==0
% %         tail = [tail(:,1)+x_offset tail(:,2)+y_offset];
% %         wholeframe = insertShape(wholeframe,'FilledCircle',[tail 3],'LineWidth',2,'Color','red');
% %         curr_tv = tail_vector(curr_idx,:).*vizscale;
% %         if sum(isnan(curr_tv))==0
% %             %             wholeframe = insertShape(wholeframe,'Line',[tail tail+curr_tv],'LineWidth',2,'Color','red');
% %         end
% %     end
    % insert bbox (added 11/10/2018)
    bbox = track.bbox(curr_idx,:);
    bbox = [bbox(1)+x_offset bbox(2)+y_offset bbox(3) bbox(4)];
%     wholeframe = insertShape(wholeframe, 'rectangle', bbox, 'Color','magenta','LineWidth',2);
    
    % insert splines (added 11/09/2018)
    curr_spline = track.spline{curr_idx}';
    if ~isempty(curr_spline) && sum(sum(isnan(curr_spline)))==0
        curr_spline(1,:) = curr_spline(1,:)+x_offset;
        curr_spline(2,:) = curr_spline(2,:)+y_offset;
        spline_out = zeros(1,2*numel(curr_spline(1,:)));
        spline_out(1:2:end-1)=curr_spline(1,:); %interleave x and y points
        spline_out(2:2:end)=curr_spline(2,:);
        wholeframe = insertShape(wholeframe,'Line',spline_out,'LineWidth',2,'Color','cyan');
        %insert head and tail based on the spline
        head = curr_spline(:,1)'; 
        wholeframe = insertShape(wholeframe,'FilledCircle',[head 3],'LineWidth',2,'Color','green');
        tail = curr_spline(:,end)'; 
        wholeframe = insertShape(wholeframe,'FilledCircle',[tail 3],'LineWidth',2,'Color','red');
    else
        disp('NO SPLINE');
    end
    
    %display speed
    speed = [num2str(round(speed_smooth(curr_idx),2)) ' mm/sec'];
    wholeframe = insertText(wholeframe, [900 10], speed, 'BoxOpacity', 1, 'FontSize', 45);
    
    %display whether animal is moving forward or backward
    curr_for = moving_forward(curr_idx);
    if curr_for
        wholeframe = insertText(wholeframe,[900 100],'FORWARD','FontSize',45,'BoxColor','green','BoxOpacity',1);
    end
    curr_back = moving_backward(curr_idx);
    if curr_back
        wholeframe = insertText(wholeframe,[900 200],'BACKWARD','FontSize',45,'BoxColor','red','BoxOpacity',1);
    end
    
    %display head grayscale value
    gs = [num2str(round(track.head_gs(curr_idx),3)) ' (gs)'];
    %     wholeframe = insertText(wholeframe, [400 10], gs,'BoxColor','blue', 'TextColor','white','BoxOpacity', 1, 'FontSize', 45);
    
    if behavFLAG
        %display whether the animal is doing a head poke
        curr_poke = head_pokes(curr_idx);
        if curr_poke
            wholeframe = insertText(wholeframe,[600 100],'POKE','FontSize',45,'BoxColor','yellow','BoxOpacity',1);
            wait = 1;
        end
        
        if ~isempty(head_pokes_rev)
            curr_poke = head_pokes_rev(curr_idx);
            if curr_poke
                wholeframe = insertText(wholeframe,[600 100],'POKE+REV','FontSize',45,'BoxColor','yellow','BoxOpacity',1);
                wait = 1;
            end
        end
        
        %display whether the animal is doing a lawn leaving event or a lawn
        %entry event
        curr_entry = lawn_entries(curr_idx);
        curr_exit = lawn_exits(curr_idx);
        if curr_entry
            wholeframe = insertText(wholeframe,[500 100],'ENTRY!','FontSize',45,'BoxColor','green','BoxOpacity',1);
            wait = 1;
        end
        if curr_exit
            wholeframe = insertText(wholeframe,[500 100],'EXIT!','FontSize',45,'BoxColor','red','BoxOpacity',1);
            wait = 1;
        end
        
        rd = roam_dwell(curr_idx);
        if rd == 1
            wholeframe = insertText(wholeframe,[200 100],'DWELL','FontSize',45,'BoxColor','white','BoxOpacity',1);
        elseif rd == 2
            wholeframe = insertText(wholeframe,[200 100],'ROAM','FontSize',45,'BoxColor','white','BoxOpacity',1);
        else
        end
        
        ehdist = track.ev_ho_dist(curr_idx);
        wholeframe = insertText(wholeframe, [size(wholeframe,2)-300 size(wholeframe,1)-200], ['EHDIST = ' num2str(round(ehdist,1))],'FontSize',30,'BoxColor','white','BoxOpacity',1);
        if ehdist < 45
            wholeframe = insertShape(wholeframe,'FilledRectangle',[size(wholeframe,2)-100 size(wholeframe,1)-100 100 100],'Color','green');
        end
    end
    
    if os_wait
        pause(2)
    end
    
    whole_player.step(wholeframe);
%         pause();
            pause(0.05);
    
    if wait %slows down playback so you can see lawn leaving and entering
        pause(0.5);
        %         pause(2);
    end
    if os_wait
        pause(2);
    end
    
end

end

