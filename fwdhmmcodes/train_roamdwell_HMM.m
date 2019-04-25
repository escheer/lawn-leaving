function [estTR, estE] = train_roamdwell_HMM(TRACKS,binSize,ratio_cutoff,trans,emis,x_offset)
%this function trains the HMM on a set of tracks and spits out the
%estimated transition and emission probabilities.

% 1. go through all TRACKS, extract onfood bits and feed these to the
% training algorithm
finalTracks = struct();

for trk = 1:length(TRACKS)
    fa = TRACKS(trk).framesActive'; %make sure to only use frames after 20 minutes in
    onfood = TRACKS(trk).centroidinlawn;
    okframes = onfood&(fa>3600); %this ensures that track is at least 20 minutes in to the video
    if sum(okframes)<binSize*2
        continue;
    end
    %find intervals when animal is fully in lawn
    onfood_ints = get_intervals( okframes, 0 );
    
    idx = 1:size(onfood_ints,1);
    ft_idx = idx+size(finalTracks,1)-1; %add this offset to keep adding more tracks in
    for i = 1:length(idx)
        finalTracks(ft_idx(i)).speed = TRACKS(trk).speed(onfood_ints(idx(i),1):onfood_ints(idx(i),2));
        finalTracks(ft_idx(i)).angspeed = TRACKS(trk).angspeed(onfood_ints(idx(i),1):onfood_ints(idx(i),2));
        finalTracks(ft_idx(i)).frames = onfood_ints(idx(i),1):onfood_ints(idx(i),2);
    end
end

statemap = getStateAuto2(finalTracks,binSize,ratio_cutoff,x_offset);

for (i=1:length(statemap)) % clean up first couple entries
    if ~isempty(statemap(i).state)
        if length(statemap(i).state) == 1 %don't include these
        else
            tmp = statemap(i).state(2:(length(statemap(i).state))); %why do 2->end? first measurement is noisy
%             tmp = statemap(i).state; %just use the whole thing 02-05-2018
            newseq(i).states = tmp;
        end

    end
end
if ~exist('newseq','var') %if we never encountered a track bit longer than a single bin, return NaNs as before
    estTR = NaN;
    estE = NaN;
    return;
end

seqs = struct2cell(newseq);
[estTR,estE] = hmmtrain(seqs,trans,emis,'MAXITERATIONS',1e4,'VERBOSE',true);
