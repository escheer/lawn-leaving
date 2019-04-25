function playframes_lawnexits( summary_struct, vidindex, indices, behavFLAG )
%PLAYFRAMES_LAWNEXITS.M This little function is a wrapper for videoplayer to watch
%defined set of frames using the summary_struct generated by
%merge_and_align_LL.m

coherence_thresh = 90; %default parameters
speed_thresh = 0.02;

% get the video you want to watch
track = summary_struct(vidindex);
if ischar(indices)
    if strcmp(indices,'all')
        indices = 1:size(track.FRAMESACTIVE,2);
    else
        error('if string is specified, must be "all" ');
    end
end

whole_player = vision.DeployableVideoPlayer('Size','Custom','CustomSize',[700 700]);

% change directory to the folder with the movies associated with this track
sourcefile = track.sourcefile;
vidfolder = uigetdir(pwd,['Find original: ' sourcefile]);
cd(vidfolder);

for idx = 1:size(track.LL_INTS,1)
    FA = track.FRAMESACTIVE(idx,:);
    CENTROID = [track.CENTX(idx,:)' track.CENTY(idx,:)'];
    HEAD = [track.HEADX(idx,:)' track.HEADY(idx,:)'];
    TAIL = [track.TAILX(idx,:)' track.TAILY(idx,:)'];
    SPEED = track.SPEED(idx,:);
    [head_vector, cent_vector, tail_vector, speed_smooth, moving_forward, moving_backward] = getforwardreverse( CENTROID, HEAD, TAIL, SPEED, coherence_thresh, speed_thresh );
    
    if behavFLAG
        head_pokes_all = track.HEADPOKES_ALL(idx,:);
        head_pokes_rev = track.HEADPOKES_REVERSAL(idx,:);
        lawn_entries = track.LAWNENTRIES(idx,:);
        lawn_exits = track.LAWNEXITS(idx,:);
        roam_dwell = track.ROAMINGDWELLINGHMM(idx,:);
    end
    
    VIDEONAME = track.VIDEONAME(idx,:);
    nonemptyvididx = find(~cellfun(@isempty,VIDEONAME));
    currvid = VIDEONAME{nonemptyvididx(1)}; %open the video corresponding to the first non-empty frame
    v = VideoReader(currvid);
    
    BGVIDINDEX = track.BGVIDINDEX(idx,:); %just for this exit
    
    last_x_offset = track.bgstruct(BGVIDINDEX(nonemptyvididx(1))).region_rounded(1);
    last_y_offset = track.bgstruct(BGVIDINDEX(nonemptyvididx(1))).region_rounded(2);
    
    for i = 1:length(indices)
        os_wait = 0;
        wait = 0;
        curr_idx = indices(i);
        disp(curr_idx);
        
%         if curr_idx == 167
%             disp('debug');
%         end
        
        %02/21/2018 there appears to be an off by one bug in BGVIDINDEX -- not
        %obvious where it comes from but just look to the next one for now as
        %long as not at the end.
        
        if curr_idx<length(BGVIDINDEX)
            next = curr_idx+1;
        else
            next = curr_idx;
        end
        
        if isnan(FA(curr_idx)) || isnan(FA(next)) %if there is no data on this frame, skip it
            disp('frame not active -- SKIP!')
            continue;
        end
        
        videoname = track.VIDEONAME{idx,curr_idx};
        videoframenum = track.VIDEOFRAME(idx,curr_idx);
        if ~strcmp(currvid,videoname)
            clear v;
            v = VideoReader(videoname); %#ok<TNMLP>
        end
        wholeframe = read(v,videoframenum);
        wholeframe = insertText(wholeframe, [10 10], curr_idx, 'BoxOpacity', 1, 'FontSize', 45);
        
        x_offset = track.bgstruct(BGVIDINDEX(next)).region_rounded(1);
        y_offset = track.bgstruct(BGVIDINDEX(next)).region_rounded(2);
        
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
        eh = track.bgstruct(BGVIDINDEX(next)).ev_ho;
        ev_ho = zeros(1,2*length(eh));
        ev_ho(1:2:length(ev_ho))=eh(:,1);
        ev_ho(2:2:length(ev_ho))=eh(:,2);
        wholeframe = insertShape(wholeframe,'Line',ev_ho,'LineWidth',2,'Color','w');
        
        %plot centroid and speed
        centroid = CENTROID(curr_idx,:);
        centroid = [centroid(:,1)+x_offset centroid(:,2)+y_offset];
        wholeframe = insertShape(wholeframe,'FilledCircle',[centroid 5],'LineWidth',2,'Color','magenta');
        
        % insert head, centroid, and tail movement vectors
        vizscale = 5;
        curr_cv = cent_vector(curr_idx,:).*vizscale;
        if sum(isnan(curr_cv))==0
            %         wholeframe = insertShape(wholeframe,'Line',[centroid centroid+curr_cv],'LineWidth',2,'Color','magenta');
        end
        head = HEAD(curr_idx,:);
        if ~isempty(head) && sum(isnan(head))==0
            head = [head(:,1)+x_offset-5 head(:,2)+y_offset-5]; %-5 until this is fixed in the original tracking code
            wholeframe = insertShape(wholeframe,'FilledCircle',[head 3],'LineWidth',2,'Color','green');
            curr_hv = head_vector(curr_idx,:).*vizscale;
            if sum(isnan(curr_hv))==0
                %             wholeframe = insertShape(wholeframe,'Line',[head head+curr_hv],'LineWidth',2,'Color','green');
            end
        end
        tail = TAIL(curr_idx,:);
        if ~isempty(tail) && sum(isnan(tail))==0
            tail = [tail(:,1)+x_offset-5 tail(:,2)+y_offset-5];
            wholeframe = insertShape(wholeframe,'FilledCircle',[tail 3],'LineWidth',2,'Color','red');
            curr_tv = tail_vector(curr_idx,:).*vizscale;
            if sum(isnan(curr_tv))==0
                %             wholeframe = insertShape(wholeframe,'Line',[tail tail+curr_tv],'LineWidth',2,'Color','red');
            end
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
        HEADGS = track.HEADGS(idx,:);
        gs = [num2str(round(HEADGS(curr_idx),3)) ' (gs)'];
        %     wholeframe = insertText(wholeframe, [400 10], gs,'BoxColor','blue', 'TextColor','white','BoxOpacity', 1, 'FontSize', 45);
        
        if behavFLAG
            %display whether the animal is doing a head poke
            curr_poke = head_pokes_all(curr_idx);
            if curr_poke
                wholeframe = insertText(wholeframe,[600 100],'POKE(ALL)','FontSize',45,'BoxColor','yellow','BoxOpacity',1);
                wait = 1;
            end
            curr_poke = head_pokes_rev(curr_idx);
            if curr_poke
                wholeframe = insertText(wholeframe,[600 100],'POKE+REV','FontSize',45,'BoxColor','yellow','BoxOpacity',1);
                wait = 1;
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
        end
        
        if os_wait
            pause(2)
        end
        
        whole_player.step(wholeframe);
%         pause(0.1);
        
        if wait %slows down playback so you can see lawn leaving and entering
            pause(0.5);
            %         pause(2);
        end
        if os_wait
            pause(2);
        end
        
    end
    pause(); %wait for user input before start the next lawn exit bout
end
end
