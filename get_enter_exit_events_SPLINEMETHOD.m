function [TRACKS, EXIT_STRUCT] = get_enter_exit_events_SPLINEMETHOD( TRACKS, before_sec, after_sec , bg_struct)
%GET_ENTER_EXIT_EVENTS_SPLINEMETHOD.M This function takes in a set of tracks and
%identifies the indices of lawn exits and entries .
%   Also returns a set of durations outside the lawn on a per track basis.
%   Also returns a set of indices surrounding each exit event, which can be
%   used to extract other information corresponding to behavior at these
%   times like speed, position, etc.

FRAMES_IN_LAWN = [];%get the frames in which the animal was in the lawn expressed in overall frame indices
ENTER_FRAMES = [];
ENTER_TRKS = [];
EXIT_TRKS = [];
EXIT_FRAMES = [];
EXIT_INTS = [];
EXIT_INT_TRK_KEY = [];
INTS_IN_BY_TRACK = [];
INTS_OUT_BY_TRACK = [];
INTS_IN_GLOBAL = [];
INTS_OUT_GLOBAL = [];
IN_INT_TRK_KEY = [];
OUT_INT_TRK_KEY = [];

for trk = 1:length(TRACKS)
%     disp(trk);
%     if trk==16
%         disp('debug');
%     end
    track = TRACKS(trk);
    bgvidindex = track.bgvidindex;
    
    lawn_entries = zeros(track.age,1);
    lawn_exits = zeros(track.age,1);
    
    % new method using entire spline 01/15/19
    spline = track.spline;
    splineDoesNotExist = cellfun(@(x) isequal(size(x),[1 1]),spline,'UniformOutput',true);
    spline(splineDoesNotExist) = deal({[NaN NaN]});
    
    splinePtsInLawn = cell(size(spline));
    for i = 1:track.age
        ev_ho = bg_struct(bgvidindex(i)).ev_ho_crp_rel; %get the event horizon that corresponds to this frame of the video
        currspline = spline{i};
        splinePtsInLawn(i) = {inpolygon(currspline(:,1),currspline(:,2),ev_ho(:,1),ev_ho(:,2))};
    end
    
    ALLil = cellfun(@(x) sum(x)==size(x,1),splinePtsInLawn,'UniformOutput',true);
    ALLol = cellfun(@(x) sum(x)==0,splinePtsInLawn,'UniformOutput',true);
    
    %remove data from frames when the worm is not segmented
    ALLil(splineDoesNotExist)=false;
    ALLol(splineDoesNotExist)=false;
        
    cil = track.centroidinlawn;
    FRAMES_IN_LAWN = [FRAMES_IN_LAWN; find(cil)+track.framesActive(1)-1];    

    % look for lawn entry as paired decrement of ALLol and increment of ALLil,
    %vice versa for lawn exit
    inside_to_X = find([0 diff(ALLil)']==-1)';
    X_to_inside = find([0 diff(ALLil)']==1)';
    outside_to_X = find([0 diff(ALLol)']==-1)';
    X_to_outside = find([0 diff(ALLol)']==1)';
    
    % by convention, we refer to leaving events as +1 and entering as -1
    glom =  [[inside_to_X ones(size(inside_to_X))];...
        [X_to_outside ones(size(X_to_outside))];
        [outside_to_X -1*ones(size(outside_to_X))];
        [X_to_inside -1*ones(size(X_to_inside))]];
    
    glom = sortrows(glom,1);
    consec = [1;diff(glom(:,2))]==0;
    cross_inds = [glom(consec,1) glom(consec,2)]; %look for consecutive upsteps and downsteps
    
    % GET THE INDICES OF LAWN ENTRY AND EXIT (TRACK INDEX)
    enter_ind = cross_inds(cross_inds(:,2)==-1,1);%two consecutive downsteps = entering
    exit_ind = cross_inds(cross_inds(:,2)==1,1); %two consecutive upsteps = exiting
    
    if ~isempty(enter_ind)
        ENTER_TRKS = [ENTER_TRKS; repmat(trk,length(enter_ind),1) enter_ind];
    end
    
    if ~isempty(exit_ind)
        EXIT_TRKS = [EXIT_TRKS; repmat(trk,length(exit_ind),1) exit_ind];
    end
    
    %ADD THESE EVENTS BACK TO THE TRACKS THEMSELVES
    lawn_entries(enter_ind) = 1;
    TRACKS(trk).lawn_entries = lawn_entries;
    
    lawn_exits(exit_ind) = 1;
    TRACKS(trk).lawn_exits = lawn_exits;
    
    % GET INTERVALS SURROUNDING LAWN LEAVING EVENTS, WHICH CAN BE USED TO
    % EXTRACT OTHER INFORMATION
    % EXTRACT AND ALIGN LAWN LEAVING EVENTS
    after = (after_sec*3)-1; %45 seconds before and 15 seconds after
    before = (before_sec*3);
    ALIGNED_EXIT_IND = before+1; %the exit_ind sits at index 136
    chopped_exits = [exit_ind-before exit_ind+after];
    
    %abbreviate these intervals at the ends of tracks and whenever another
    %entry or exit event is encountered -- pad with NaNs in these indices.
    if ~isempty(chopped_exits)
        for i = 1:size(chopped_exits,1)
            ext = exit_ind(i);
            int = chopped_exits(i,1):chopped_exits(i,2);
            inrange = int(int>0 & int<=track.age);
            abbrev = inrange;
            badboys = sort(intersect(union(setdiff(exit_ind,ext),enter_ind),inrange)); %check if these intervals run into subsequent lawn leaving or entering events or the end or beginning of the track
            if ~isempty(badboys)
                upperlim = min(badboys(badboys>ext))-1;
                lowerlim = max(badboys(badboys<ext))+1;
                if ~isempty(upperlim)
                    abbrev = abbrev(1):upperlim;
                end
                if ~isempty(lowerlim)
                    abbrev = lowerlim:abbrev(end);
                end
            end
            [~,ia] = setdiff(int,abbrev); %replace out of bounds indices from int with NaNs, then add it
            int(ia)=NaN;
            EXIT_INTS = [EXIT_INTS; int];
            EXIT_INT_TRK_KEY = [EXIT_INT_TRK_KEY; trk];
        end
    end
    ENTER_FRAMES = [ENTER_FRAMES; enter_ind+track.framesActive(1)-1];
    EXIT_FRAMES = [EXIT_FRAMES; exit_ind+track.framesActive(1)-1]; %expressed in terms of total frames in the video
    
    %FIGURE OUT THE INTERVALS IN AND OUT
    if ~(isempty(enter_ind) && isempty(exit_ind)) %if there is at least a single enter or exit event in the track, proceed
        tmp_entry = enter_ind; tmp_exit = exit_ind;
        if min(enter_ind)<min(exit_ind) %worm began outside lawn, then entered, so 1st frame acts as an exit
            tmp_exit = union(1,tmp_exit);
            started_out = 1;
        elseif min(exit_ind)<min(enter_ind) %worm began inside lawn, then exited, so 1st frame acts as an entry
            tmp_entry = union(1,tmp_entry);
            started_out = 0;
        elseif isempty(enter_ind)
            tmp_entry = union(1,tmp_entry);
            started_out = 0;
        elseif isempty(exit_ind)
            tmp_exit = union(1,tmp_exit);
            started_out = 1;
        else
            error('no other option!');
        end
        if max(enter_ind)<max(exit_ind) %worm began outside lawn, then entered, so 1st frame acts as an exit
            tmp_entry = union(tmp_entry, track.age);
        elseif max(exit_ind)<max(enter_ind) %worm began inside lawn, then exited, so 1st frame acts as an entry
            tmp_exit = union(tmp_exit, track.age);
        elseif isempty(enter_ind)
            tmp_entry = union(tmp_entry, track.age);
        elseif isempty(exit_ind)
            tmp_exit = union(tmp_exit, track.age);
        else
            error('no other option!');
        end
        
        total_events = length(tmp_entry)+length(tmp_exit);
        t = zeros(1,total_events);
        twoways = logical(toeplitz(mod(1:total_events,2),mod(1:2,2))); %column 1 is alternating starting with 1, column 2 with 0, 1 means you're OUT
        if started_out %started out so index 1 means youre out of the lawn
            t(1:2:total_events) = tmp_exit;
            t(2:2:total_events) = tmp_entry;
            in_or_out = twoways(:,1);
        else %started in so index 1 means youre in the lawn
            t(1:2:total_events) = tmp_entry;
            t(2:2:total_events) = tmp_exit;
            in_or_out = twoways(:,2);
        end
        
        for i = 1:length(t)-1
            first = t(i); last = t(i+1)-1;
            if i == length(t)-1
                last = t(i+1);
            end
            if in_or_out(i)
                INTS_OUT_BY_TRACK = [INTS_OUT_BY_TRACK; first last];
                OUT_INT_TRK_KEY = [OUT_INT_TRK_KEY; trk];
                INTS_OUT_GLOBAL = [INTS_OUT_GLOBAL; first+track.framesActive(1)-1 last+track.framesActive(1)-1];
            else
                INTS_IN_BY_TRACK = [INTS_IN_BY_TRACK; trk first last];
                IN_INT_TRK_KEY = [IN_INT_TRK_KEY; trk];
                INTS_IN_GLOBAL = [INTS_IN_GLOBAL; first+track.framesActive(1)-1 last+track.framesActive(1)-1];
            end
        end
    else
        first = 1; last = track.age;
        if sum(track.fullyinlawn) == 0 %this track neither enters nor exits; it is OUT the entire time
            INTS_OUT_BY_TRACK = [INTS_OUT_BY_TRACK; first last];
            OUT_INT_TRK_KEY = [OUT_INT_TRK_KEY; trk];
            INTS_OUT_GLOBAL = [INTS_OUT_GLOBAL; first+track.framesActive(1)-1 last+track.framesActive(1)-1];
        else
            INTS_IN_BY_TRACK = [INTS_IN_BY_TRACK; trk first last];
            IN_INT_TRK_KEY = [IN_INT_TRK_KEY; trk];
            INTS_IN_GLOBAL = [INTS_IN_GLOBAL; first+track.framesActive(1)-1 last+track.framesActive(1)-1];
        end
    end
end

EXIT_STRUCT.FRAMES_IN_LAWN = FRAMES_IN_LAWN;
EXIT_STRUCT.ENTER_FRAMES = ENTER_FRAMES;
EXIT_STRUCT.ENTER_TRKS = ENTER_TRKS;
EXIT_STRUCT.EXIT_FRAMES = EXIT_FRAMES;
EXIT_STRUCT.EXIT_TRKS = EXIT_TRKS;
EXIT_STRUCT.EXIT_INTS = EXIT_INTS;
EXIT_STRUCT.ALIGNED_EXIT_IND = ALIGNED_EXIT_IND;
EXIT_STRUCT.EXIT_TRK_KEY = EXIT_INT_TRK_KEY;
EXIT_STRUCT.INTS_IN_BY_TRACK = INTS_IN_BY_TRACK;
EXIT_STRUCT.INTS_OUT_BY_TRACK = INTS_OUT_BY_TRACK;
EXIT_STRUCT.INTS_IN_GLOBAL = INTS_IN_GLOBAL;
EXIT_STRUCT.INTS_OUT_GLOBAL = INTS_OUT_GLOBAL;
EXIT_STRUCT.IN_INT_TRK_KEY = IN_INT_TRK_KEY;
EXIT_STRUCT.OUT_INT_TRK_KEY = OUT_INT_TRK_KEY;

