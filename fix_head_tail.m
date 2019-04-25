function TRACKS = fix_head_tail( TRACKS, bg_struct )
%FIX_HEAD_TAIL.m This function fixes errors in head/tail annotation after
%analysis. Re-runs smoothing function to form head_smooth and tail_smooth.
%Re-orders head and tail grayscale values accordingly. Re-runs getforward
%Will eventually be superseded by fixing these errors in
%head_and_tail_kalman.m But will remain an essential check before
%proceeding on to other analyses that derive from information about
%head/tail annotation.
for trk = 1:length(TRACKS)
    track = TRACKS(trk);
    head = track.head;
    tail = track.tail;
    head_gs = track.head_gs;
    tail_gs = track.tail_gs;
    
    %look for discontinuities in smooth head/tail annotation -- these often
    %arise from poor segmentation due to a piece of dirt or a laid egg.
    head_d = diag(pdist2(head(1:end-1,:),head(2:end,:)));% distance from head at frame x to frame x+1
    tail_d = diag(pdist2(tail(1:end-1,:),tail(2:end,:)));
    %re-examine any intervals bookended by a jump of more than 20 pixels to see
    %if it is heads or tails
    thresh = 20;
    a = find(head_d>thresh | tail_d>thresh);
    
    if ~isempty(a)
        split_starts = a((diff(a)~=1)');
        split_starts = unique([split_starts;a(end)]);%add back last zero member, sometimes lost
        a_diff = (diff(a)==1)';
        f = find([false,a_diff]~=[a_diff,false]);
        g = find(f(2:2:end)-f(1:2:end-1));
        split_ends = a(f(2*g-1));
        splits = unique(sort([1;split_starts;split_ends;length(head_d)])); %these are all breakpoints
        check_int = [splits(1:end-1)+1 splits(2:end)-1];
        check_int(1) = 1; check_int(end) = length(head_d); %make sure the first index is 1 and the last index is the last number in the array
        check_int(check_int(:,1)-check_int(:,2)>0,:) = []; %get rid of impossible intervals
        % and intervals to keep, those containing nonzero, nonnan entries
    else %no need to split this track further
        check_int = [1 track.totalVisibleCount];
    end
    
    %for any of the intervals less than 10 frames long, these will be deleted
    %as they are too ambiguous based on movement alone
    int_dur = check_int(:,2)-check_int(:,1);
    delete_int = check_int(int_dur<10,:);
    check_int(int_dur<10,:) = []; %remove these from consideration
    
    %check which is really the head among check_int
    for k = 1:size(check_int,1)
        idx = check_int(k,1):check_int(k,2);
        washead = head(idx,:);
        washead_gs = head_gs(idx);
        wastail = tail(idx,:);
        wastail_gs = tail_gs(idx);
        
        a_H = squareform(pdist(washead));
        b_H = a_H(:,2:end);
        c_H = diag(b_H);
        d_H = nansum(c_H);
        
        a_T = squareform(pdist(wastail));
        b_T = a_T(:,2:end);
        c_T = diag(b_T);
        d_T = nansum(c_T);
        
        if d_H-d_T>=5 % head remains the head, tail remains tail
            head(idx,:) = washead;
            tail(idx,:) = wastail;
            head_gs(idx) = washead_gs;
            tail_gs(idx) = wastail_gs;
        elseif d_T-d_H>=5 % head becomes tail, tail becomes head
            head(idx,:) = wastail;
            tail(idx,:) = washead;
            head_gs(idx) = wastail_gs;
            tail_gs(idx) = washead_gs;
            disp(['TRACK '  num2str(trk) ' frame ' num2str(idx(1)) 'to' num2str(idx(end)) ': HEAD AND TAIL CORRECTED!']);
        else
            head(idx,:) = NaN(size(idx,2),2);
            tail(idx,:) = NaN(size(idx,2),2);
            head_gs(idx) = NaN(size(idx,2),1);
            tail_gs(idx) = NaN(size(idx,2),1);
            warning(['frame ' num2str(idx(1)) 'to' num2str(idx(end)) ': cannot discern H v T by displacement in this interval!']);
        end
    end
    
    %remove aforementioned ambiguous intervals
    for k = 1:size(delete_int,1)
        idx = delete_int(k,1):delete_int(k,2);
        head(idx,:) = NaN(size(idx,2),2);
        tail(idx,:) = NaN(size(idx,2),2);
        
        head_gs(idx,:) = NaN(size(idx,2),1);
        tail_gs(idx,:) = NaN(size(idx,2),1);
    end
    %re-calculate smoothed position, this time window is 3
    head_smooth = [movmean(head(:,1),3,'omitnan') movmean(head(:,2),3,'omitnan')];
    tail_smooth = [movmean(tail(:,1),3,'omitnan') movmean(tail(:,2),3,'omitnan')];
    
    %re-calculate forward and backward intervals based on corrected head
    %and tail annotation.
    [~, ~, ~, speed_smooth, forward, reverse] = getforwardreverse( track.centroid,head,tail,track.speed, 90, 0.02 );
        
    %update TRACKS
    TRACKS(trk).head = head;
    TRACKS(trk).head_smooth = head_smooth;
    TRACKS(trk).tail = tail;
    TRACKS(trk).tail_smooth = tail_smooth;
    TRACKS(trk).head_gs = head_gs;
    TRACKS(trk).tail_gs = tail_gs;
    TRACKS(trk).speed_smooth = speed_smooth;
    TRACKS(trk).forward = forward;
    TRACKS(trk).reverse = reverse;
    
    %re-calculate when centroid, head, and tail are in the lawn, update
    %TRACKS (this is done in this order because the function
    %countBlobsInOut reads values out of the track as it is saved
    [headinlawn, centroidinlawn, tailinlawn, fullyinlawn] = countBlobsInOut(TRACKS(trk), bg_struct);
    TRACKS(trk).headinlawn = headinlawn;
    TRACKS(trk).centroidinlawn = centroidinlawn;
    TRACKS(trk).tailinlawn = tailinlawn;
    TRACKS(trk).fullyinlawn = fullyinlawn;
    
end

