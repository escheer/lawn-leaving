function tracks = head_and_tail( tracks )
%head_and_tail.m This function simply tracks the ends of the fitted spline
%based on which one is closer through time and returns the two tracks. In
%cases where the head and tail are both very close (i.e. a reorientation),
%it breaks the tracks. Then using the color metric, it re-sorts the tracks
%into head and tail chunks.

side1 = {};
side2 = {};
subtrack_counter = 0;
for trackidx = 1:length(tracks)
    disp(['subtrack counter' num2str(trackidx)]);
    %     if trackidx==2
    %         disp('debug');
    %     end
    subtrack_counter = subtrack_counter+1; %increment cell counter for new tracks
    side1{subtrack_counter} = [];
    side2{subtrack_counter} = [];
    last_NaN = 1;
    spline = tracks(trackidx).spline;
    split_track = 0;
    %     last_NaN = 1; %always start fresh when a new track starts
    for i = 1:length(spline)
        disp(['spline counter ' num2str(i)]);
        curr_spline = spline{i};
        if ~isnan(curr_spline)
            if last_NaN %initialize by random choice if it's the first one or if there have just been NaNs
                if i > 1
                    subtrack_counter = subtrack_counter+1; %will have just incremented this if its a new track
                end
                side1{subtrack_counter} = [curr_spline(1,:) trackidx i];
                side2{subtrack_counter} = [curr_spline(end,:) trackidx i];
                last_NaN = 0;
            else
                last_side1 = side1{subtrack_counter}(end,:);
                last_side2 = side2{subtrack_counter}(end,:);
                endpt1 = curr_spline(1,:);
                endpt2 = curr_spline(end,:);
                d = pdist2([last_side1(1:2);last_side2(1:2)],[endpt1;endpt2]);
                [d_1,idx1]=min(d(1,:));
                [d_2,idx2]=min(d(2,:));
                if idx1 == idx2 %both previous points are closer to one of the current points than the other
                    if d_1 < d_2
                        idx1 = 1;
                        idx2 = 2;
                    elseif d_2 > d_1
                        idx1 = 2;
                        idx2 = 1;
                    else %if d2 == d1
                        split_track = 1; %send to the split track clause
                        disp('Both previous points are closer to one of the new endpoints');
                        %                         error('idx 1 cannot be same as idx 2!');
                    end
                end
                if (~abs(d(1,1)-d(1,2))<10 || ~abs(d(2,1)-d(2,2))<10) && ~split_track
                    if idx1 == 1 && idx2 == 2
                        side1{subtrack_counter} = [side1{subtrack_counter}; [endpt1 trackidx i]];
                        side2{subtrack_counter} = [side2{subtrack_counter}; [endpt2 trackidx i]];
                    elseif idx1 == 2 && idx2 == 1
                        side1{subtrack_counter} = [side1{subtrack_counter}; [endpt2 trackidx i]];
                        side2{subtrack_counter} = [side2{subtrack_counter}; [endpt1 trackidx i]];
                    else
                        error('UNEXPECTED!');
                    end
                    last_NaN = 0;
                else %too close too call, split track
                    disp('----AMBIGUOUS: SPLIT TRACK----');
                    side1{subtrack_counter} = [side1{subtrack_counter}; [NaN NaN trackidx i]];
                    side2{subtrack_counter} = [side2{subtrack_counter}; [NaN NaN trackidx i]];
                    warning(['frame: ' num2str(i) ' h v t too close to call, split track! reorientation???']);
                    last_NaN = 1;
                    split_track = 0;
                end
            end
        else
            disp('spline is NaN; put NaN and split track!');
            side1{subtrack_counter} = [side1{subtrack_counter}; [NaN NaN trackidx i]];
            side2{subtrack_counter} = [side2{subtrack_counter}; [NaN NaN trackidx i]];
            last_NaN = 1;
        end
    end
end
%now for each pair of subtracks, decide which one is head or tail by
%brightness test, stitch heads and tails back together
disp('----------NOW STITCH BACK TOGETHER---------');
head = [];
tail = [];
last_track = 0;
debug = 0;
for j = 1:length(side1)
    debug = 0;
    disp(['side1 loop ' num2str(j)]);
%     if j==14
%         disp('debug');
%         debug = 1;
%     end
    s1 = side1{j};
    s2 = side2{j};
    curr_track = s1(1,3);
    if curr_track == last_track+1 %re-initialize for each new track
        head = [];
        tail = [];
        last_track = curr_track;
    end
    s1_inds = NaN(size(s1,1),1);
    s2_inds = NaN(size(s2,1),1);
    for k = 1:size(s1,1)
%         if debug
%             disp(k)
% %             if k == 29
% %                 disp('debug');
% %             end
%         end
        curr_s1 = s1(k,:);
        curr_s2 = s2(k,:);
        if ~isnan(curr_s1(1)) || ~isnan(curr_s2(1))
            wormcrop_bin = tracks(curr_s1(3)).cropworm{curr_s1(4)};
            wormcrop_orig = tracks(curr_s1(3)).cropworm_orig{curr_s1(4)};
            bbx = tracks(curr_s1(3)).bbox(curr_s1(4),:);
            x_offset = bbx(1);
            y_offset = bbx(2);
            endpt_1 = [curr_s1(1)-x_offset curr_s1(2)-y_offset];
            endpt_2 = [curr_s2(1)-x_offset curr_s2(2)-y_offset];
            [e1_call, e2_call] = get_worm_head(endpt_1,endpt_2,wormcrop_bin,wormcrop_orig);
            s1_inds(k) = e1_call;
            s2_inds(k) = e2_call;
        else
            s1_inds(k) = NaN;
            s2_inds(k) = NaN;
        end
    end

%     s1_call = mode(s1_inds);
%     s2_call = mode(s2_inds);
    if sum(isnan(s1_inds))==length(s1_inds) || sum(isnan(s2_inds))==length(s2_inds)
        disp('this section is filled with NaNs!');
        head = [head; NaN(size(s1))];
        tail = [tail; NaN(size(s2))];
    elseif sum(s1_inds==1)/sum(~isnan(s1_inds))>0.70 && sum(s2_inds==2)/sum(~isnan(s2_inds))>0.70
        disp('s1 is the head, s2 is the tail');
        head = [head; s1];
        tail = [tail; s2];
    elseif sum(s1_inds==2)/sum(~isnan(s1_inds))>0.70 && sum(s2_inds==1)/sum(~isnan(s2_inds))>0.70
        disp('s2 is the head, s1 is the tail');
        head = [head; s2];
        tail = [tail; s1];
    else
        disp('uncertain head and tail annotation based on 70% cutoff based on brightness test');
        head = [head; NaN(size(s1))];
        tail = [tail; NaN(size(s2))];
%         error('IDK HOW TO SAY WHICH ONE THE END IS?');
    end
    head_so_far = tracks(curr_track).head;
    if sum(isnan(head_so_far))==1
        head_so_far = [];
    end
    tracks(curr_track).head = head(:,1:2);
    tail_so_far = tracks(curr_track).tail;
    if sum(isnan(tail_so_far))==1
        tail_so_far = [];
    end
    tracks(curr_track).tail = tail(:,1:2);
    
end


end