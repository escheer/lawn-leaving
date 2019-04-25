function [expNewSeq_out, expStates_out, estTR, estE, actualratio_out] = getHMMStates2(track,binSize,onfood,ratio,trans,emis)
%function [estTR estE] = getHMMStates(finalTracks,binSize, trans,emis)

% if sum(onfood)<binSize*2 %only analyze for tracks with onfood at least 2*binSize frames
%     expNewSeq_out = NaN(size(onfood));
%     expStates_out = NaN(size(onfood));
%     estTR = NaN;
%     estE = NaN;
%     actualratio_out = NaN(size(onfood));
%     return;
% end
% 
% %chop up trackData into bits that are on food, and do the analysis on those
% % fill the rest with NaNs
% a = find(onfood==0);
% onfood_ints = [];
% if ~isempty(a)
%     split_starts = a((diff(a)~=1)');
%     split_starts = unique([split_starts;a(end)]);%add back last zero member, sometimes lost
%     a_diff = (diff(a)==1)';
%     f = find([false,a_diff]~=[a_diff,false]);
%     g = find(f(2:2:end)-f(1:2:end-1));
%     split_ends = a(f(2*g-1));
%     splits = unique(sort([1;split_starts;split_ends;length(onfood)])); %these are all breakpoints
%     split_int = [splits(1:end-1)+1 splits(2:end)-1];
%     split_int(1) = 1; split_int(end) = length(onfood); %make sure the first index is 1 and the last index is the last number in the array
%     split_int(split_int(:,1)-split_int(:,2)>0,:) = []; %get rid of impossible intervals
%     % and intervals to keep, those containing nonzero, nonnan entries
%     onfood_ints = [];
%     for n = 1:size(split_int,1)
%         c = onfood(split_int(n,1):split_int(n,2));
%         if sum(c==0)~=length(c)
%             onfood_ints = [onfood_ints ; split_int(n,:)];
%         end
%     end
% else %no need to split this track further
%     onfood_ints = [1 length(onfood)];
% end
% finalTracks = struct();
% for idx = 1:size(onfood_ints,1)
%     finalTracks(idx).speed = track.speed(onfood_ints(idx,1):onfood_ints(idx,2));
%     finalTracks(idx).angspeed = track.angspeed(onfood_ints(idx,1):onfood_ints(idx,2));
%     finalTracks(idx).frames = onfood_ints(idx,1):onfood_ints(idx,2);
% end

%TRY WITHOUT ONFOOD STIPULATION:
finalTracks = struct();
finalTracks.speed = track.speed;
finalTracks.angspeed = track.angspeed;
finalTracks.frames = 1:track.age;

% call this for every contiguous bit of inlawn track

%what to do if the length of finalTracks is less than the bin size? FIX
%THIS 01-30-2018%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

statemap = getStateAuto(finalTracks,binSize,ratio);

statesHMM = [];
for (i=1:length(statemap)) % clean up first couple entries
    if ~isempty(statemap(i).state)
        if length(statemap(i).state) == 1 %don't include these
        else
            tmp = statemap(i).state(2:(length(statemap(i).state)));
            newseq(i).states = tmp;
        end

    end
end
if ~exist('newseq','var') %if we never encountered a track bit longer than a single bin, return NaNs as before
    expNewSeq_out = NaN(size(onfood));
    expStates_out = NaN(size(onfood));
    estTR = NaN;
    estE = NaN;
    actualratio_out = NaN(size(onfood));
    return;
end

seqs = struct2cell(newseq);

%%% FOR TWO STATE MODEL
% trans = [0.995, 0.005; 0.07, 0.93];
% emis = [0.96, 0.04; 0.07, 0.93];

%%%  FOR ONE STATE MODEL
%trans = 1;
%emis = [0.8, 0.2];


%%% FOR THREE STATE MODEL
%trans = [0.98 0.01 0.01; 0.33 0.33 0.33; 0.03 0.03 0.96];
%emis = [0.96, 0.04; 0.5 0.5; 0.07, 0.93];


%%%%TRAIN THE MODEL
[estTR,estE] = hmmtrain(seqs,trans,emis); %might not need to do this every time, just use what you find from a big N2 population

%%%%Apply Viterbi algorithm
for i = 1:length(newseq)
        statesHMM(i).states = hmmviterbi(newseq(i).states,estTR,estE);
        statesHMM(i).actualratio = statemap(i).actualratio(1:end-1);
end


%%%Create expStates structure - most probable state path by HMM
%%%Re-normalize to track length, given binSize
for i = 1:length(newseq)
    if ~isempty(statesHMM(i).states)
        speedData = finalTracks(i).speed;
        numbBins = (length(speedData))/binSize;
        numbBins = floor(numbBins);
        expStates(i).states(1:(binSize)) = statesHMM(i).states(1);
        expStates(i).actualratio(1:(binSize)) = statesHMM(i).actualratio(1);
        for (j = 2:numbBins)
            startPl = (j * binSize) - (binSize-1);
            stopPl  = (j * binSize);
            expStates(i).states(startPl:stopPl) = statesHMM(i).states(j-1);
            expStates(i).actualratio(startPl:stopPl) = statesHMM(i).actualratio(j-1);
        end
        numFrames = length(finalTracks(i).frames);
        leftOver_Frames = mod(numFrames,binSize);
        if (leftOver_Frames>0)
            currentExpStateLength = length(expStates(i).states);
            newStart = currentExpStateLength+1;
            %             newStop = newStart + (leftOver_Frames -4);
            newStop = newStart + leftOver_Frames-1;
            expStates(i).states(newStart:newStop) = statesHMM(i).states(end);
            expStates(i).actualratio(newStart:newStop) = statesHMM(i).actualratio(end);
            %         else
            %             expStates(i).states = expStates(i).states(1:(numFrames-3));
        end
    end
end

%%%Create expNewSeq structure - just the binary data from 2D scatter
% for(i= 1:length(finalTracks))
for i = 1:length(newseq)
    if ~isempty(statesHMM(i).states)
        speedData = finalTracks(i).speed;
        numbBins = (length(speedData))/binSize;
        numbBins = floor(numbBins);
        expNewSeq(i).states(1:(binSize)) = newseq(i).states(1);
        %         expNewSeq(i).actualratio(1:(binSize)) = statesHMM(i).actualratio(1);
        for (j = 2:numbBins)
            startPl = (j * binSize) - (binSize-1);
            stopPl  = (j * binSize);
            expNewSeq(i).states(startPl:stopPl) = newseq(i).states(j-1);
            %             expNewSeq(i).actualratio(startPl:stopPl) = statesHMM(i).actualratio(j-1);
        end
        numFrames = length(finalTracks(i).frames);
        leftOver_Frames = mod(numFrames,binSize);
        if (leftOver_Frames>0)
            currentExpStateLength = length(expNewSeq(i).states);
            newStart = currentExpStateLength+1;
            %             newStop = newStart + (leftOver_Frames -4);
            newStop = newStart + leftOver_Frames-1;
            expNewSeq(i).states(newStart:newStop) = newseq(i).states(end);
            %             expNewSeq(i).actualratio(newStart:newStop) = statesHMM(i).actualratio(end);
            %         else
            %             expStates(i).states = expStates(i).states(1:(numFrames-3));
        end
    end
end

expStates_out = NaN(size(onfood));
expNewSeq_out = NaN(size(onfood));
actualratio_out = NaN(size(onfood));

for i = 1:length(newseq)
    if ~isempty(expStates(i).states)
        expStates_out(finalTracks(i).frames) = expStates(i).states';
        expNewSeq_out(finalTracks(i).frames) = expNewSeq(i).states';
        actualratio_out(finalTracks(i).frames) = expStates(i).actualratio';
    end
end
