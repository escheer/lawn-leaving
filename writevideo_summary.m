function writevideo_summary( SUMMARY_STRUCT, bg_struct, indices, behavFLAG, outname )
%WRITE_TRACK_VIDEO_SUMMARY This function writes out segment of video to a file
%corresponding to the indices in a SUMMARY_STRUCT of interest. Annotates head pokes,
%lawn entry and exit on the video. displays speed, forward and backward
%locomotion.

structlen = length(SUMMARY_STRUCT.SPEED);
if ischar(indices)
    if strcmp(indices,'all')
        indices = 1:structlen;
    else
        error('if string is specified, must be "all" ');
    end
end

v_out = VideoWriter([outname '.avi'],'Motion JPEG AVI');
v_out.FrameRate = 21; %7x sped up from 3fps
v_out.Quality = 100;
open(v_out);

% if ~deployable
%     whole_player = vision.VideoPlayer('Position', [1043, 32, 816, 593]);
% else
%     whole_player = vision.DeployableVideoPlayer('Size','Custom','CustomSize',[700 700]);
% end

speed_smooth = SUMMARY_STRUCT.SPEED_SMTH;
forward = SUMMARY_STRUCT.FORWARD;
reverse = SUMMARY_STRUCT.REVERSE;

if behavFLAG
    headpoke_fwd = SUMMARY_STRUCT.HEADPOKE_FWD;
    if isempty(headpoke_fwd)
        headpoke_fwd = zeros(structlen,1);
    end
    headpoke_rev = SUMMARY_STRUCT.HEADPOKE_REV;
    if isempty(headpoke_rev)
        headpoke_rev = zeros(structlen,1);
    end
    headpoke_pause = SUMMARY_STRUCT.HEADPOKE_PAUSE;
    if isempty(headpoke_pause)
        headpoke_pause = zeros(structlen,1);
    end
    lawn_entry = SUMMARY_STRUCT.LAWN_ENTRY;
    lawn_exit = SUMMARY_STRUCT.LAWN_EXIT;
    roam_dwell = SUMMARY_STRUCT.ROAMDWELL_HMM;
end
currvid = SUMMARY_STRUCT.VIDEONAME{indices(1)};
v = VideoReader(currvid);

for i = 1:length(indices)
    wait = 0;
    curr_idx = indices(i);
    
    disp(curr_idx);
    
    videoname = SUMMARY_STRUCT.VIDEONAME{curr_idx};
    
    if isempty(videoname)
        disp('NO DATA FOR CURRENT FRAME, DISPLAYING LAST FRAME!');
        wholeframe = insertText(wholeframe, [10 10], curr_idx, 'BoxOpacity', 1, 'FontSize', 45);
        writeVideo(v_out,wholeframe);
        continue;
    end
    
    videoframenum = SUMMARY_STRUCT.VIDEOFRAME(curr_idx);
    if ~strcmp(currvid,videoname)
        clear v;
        v = VideoReader(videoname); %#ok<TNMLP>
    end
    wholeframe = read(v,videoframenum);
    wholeframe = insertText(wholeframe, [10 10], curr_idx, 'BoxOpacity', 1, 'FontSize', 45);
    
    x_offset = bg_struct(SUMMARY_STRUCT.BGVIDINDEX(curr_idx)).region_rounded(1);
    y_offset = bg_struct(SUMMARY_STRUCT.BGVIDINDEX(curr_idx)).region_rounded(2);
    
    % plot event horizon
    eh = bg_struct(SUMMARY_STRUCT.BGVIDINDEX(curr_idx)).ev_ho;
    ev_ho = zeros(1,2*length(eh));
    ev_ho(1:2:length(ev_ho))=eh(:,1);
    ev_ho(2:2:length(ev_ho))=eh(:,2);
    wholeframe = insertShape(wholeframe,'Line',ev_ho,'LineWidth',2,'Color','w');
    
    %plot centroid and speed
    centroid = SUMMARY_STRUCT.CENTROID(curr_idx,:);
    if sum(isnan(centroid))==0
        centroid = [centroid(:,1)+x_offset centroid(:,2)+y_offset];
        wholeframe = insertShape(wholeframe,'FilledCircle',[centroid 5],'LineWidth',2,'Color','magenta');
    end
    %     % insert bbox
    %     bbox = SUMMARY_STRUCT.BBOX(curr_idx,:);
    %     bbox = [bbox(1)+x_offset bbox(2)+y_offset bbox(3) bbox(4)];
    %     wholeframe = insertShape(wholeframe, 'rectangle', bbox, 'Color','magenta','LineWidth',2);
    
    % insert splines
    curr_spline_x = SUMMARY_STRUCT.SPLINE_x(curr_idx,:);
    curr_spline_y = SUMMARY_STRUCT.SPLINE_y(curr_idx,:);
    curr_spline = [curr_spline_x; curr_spline_y];
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
    disp_speed = [num2str(round(speed_smooth(curr_idx),2)) ' mm/sec'];
    wholeframe = insertText(wholeframe, [900 10], disp_speed, 'BoxOpacity', 1, 'FontSize', 45);
    
    %display whether animal is moving forward or backward
    curr_for = forward(curr_idx);
    if ~isnan(curr_for)
        if curr_for
            wholeframe = insertText(wholeframe,[900 100],'FORWARD','FontSize',45,'BoxColor','green','BoxOpacity',1);
        end
    end
    curr_back = reverse(curr_idx);
    if ~isnan(curr_back)
        if curr_back
            wholeframe = insertText(wholeframe,[900 200],'BACKWARD','FontSize',45,'BoxColor','red','BoxOpacity',1);
        end
    end
    
    %     %display head grayscale value
    %     gs = [num2str(round(SUMMARY_STRUCT.HEAD_GS(curr_idx),3)) ' (gs)'];
    %     wholeframe = insertText(wholeframe, [400 10], gs,'BoxColor','blue', 'TextColor','white','BoxOpacity', 1, 'FontSize', 45);
    
    if behavFLAG
        %display whether the animal is doing a head poke
        curr_poke = headpoke_fwd(curr_idx);
        if ~isnan(curr_poke)
            if curr_poke
                disp('HP+FWD!');
                wholeframe = insertText(wholeframe,[600 100],'POKE+FWD','FontSize',45,'BoxColor','yellow','BoxOpacity',1);
                wait = 1;
            end
        end
        
        curr_poke = headpoke_rev(curr_idx);
        if ~isnan(curr_poke)
            if curr_poke
                disp('HPREV!');
                wholeframe = insertText(wholeframe,[600 100],'POKE+REV','FontSize',45,'BoxColor','yellow','BoxOpacity',1);
                wait = 1;
            end
        end
        curr_poke = headpoke_pause(curr_idx);
        if ~isnan(curr_poke)
            if curr_poke
                disp('HP-PAUSE!');
                wholeframe = insertText(wholeframe,[600 100],'POKE+PAUSE','FontSize',45,'BoxColor','yellow','BoxOpacity',1);
                wait = 1;
            end
        end
        %display whether the animal is doing a lawn leaving event or a lawn
        %entry event
        curr_entry = lawn_entry(curr_idx);
        curr_exit = lawn_exit(curr_idx);
        if ~isnan(curr_entry)
            if curr_entry
                wholeframe = insertText(wholeframe,[500 100],'ENTRY!','FontSize',45,'BoxColor','green','BoxOpacity',1);
                wait = 1;
            end
        end
        if ~isnan(curr_exit)
            if curr_exit
                wholeframe = insertText(wholeframe,[500 100],'EXIT!','FontSize',45,'BoxColor','red','BoxOpacity',1);
                wait = 1;
            end
        end
        
        rd = roam_dwell(curr_idx);
        if ~isnan(rd)
            if rd == 1
                wholeframe = insertText(wholeframe,[200 100],'DWELL','FontSize',45,'BoxColor','white','BoxOpacity',1);
            elseif rd == 2
                wholeframe = insertText(wholeframe,[200 100],'ROAM','FontSize',45,'BoxColor','white','BoxOpacity',1);
            else
            end
        end
        
        ehdist = SUMMARY_STRUCT.EV_HO_DIST(curr_idx);
        if ~isnan(ehdist)
            wholeframe = insertText(wholeframe, [size(wholeframe,2)-300 size(wholeframe,1)-200], ['EHDIST = ' num2str(round(ehdist,1))],'FontSize',30,'BoxColor','white','BoxOpacity',1);
            if ehdist < 45
                wholeframe = insertShape(wholeframe,'FilledRectangle',[size(wholeframe,2)-100 size(wholeframe,1)-100 100 100],'Color','green');
            end
        end
    end
    
    
    if wait %slows down playback so you can see lawn leaving and entering
        counter = 0;
        while counter<21
            writeVideo(v_out,wholeframe);
            counter = counter+1;
        end
    else
        writeVideo(v_out,wholeframe);
    end
    
    
end
close(v_out);

end

