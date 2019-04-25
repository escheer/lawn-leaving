function cent_lookback = getprevdist2(tracks,tracksThatLeft,index,curr_cent,numlookback,frame)
%GETPREVDIST.M Given a track, get the distance of the current position
%from the previous locations in order to figure out how far back to look
%for a frame without the worm for grayscale value extraction.

mindist = 15;
allTracks = [tracksThatLeft,tracks];
range = index-numlookback:index;
%figure out which frames in which tracks this is made up of
cent = NaN(size(range,2),2);
for i = 1:length(allTracks)
    fa = allTracks(i).framesActive;
    [~,IA,IB] = intersect(fa,range);
    cent(IB,:) = allTracks(i).centroid(IA,:);
end

if sum(isnan(cent(:,1)))>0.5*size(cent,1) %if a track just started or if the lookback is mostly frames when there were no tracks, just go back as far as possible
    cent_lookback = 360;
    return;
end

cent_dist = diag(pdist2(cent,repmat(curr_cent,length(cent),1)));
%find peaks
[c_pks,c_locs] = findpeaks(cent_dist,'MinPeakHeight',mindist,'MinPeakProminence',5);
c_locs_rel = c_locs+(index-numlookback);

if isempty(c_pks) %either all the peaks are too low or its just smoothly decreasing i.e. when the worm is moving fast
    [~,cent_closest] = nanmax(abs(cent_dist-(mindist+5))<=1);
    if cent_dist(cent_closest)<mindist %if the best you can do is less than mindist+5, just take the max distance
        [~, cent_closest] = nanmax(cent_dist);
    end
    cent_closest_rel = cent_closest+(index-numlookback);
    cent_lookback = index-cent_closest_rel;
%     min_cent_dist = cent_dist(cent_closest);
%     c_locs = cent_closest;
%     c_pks = min_cent_dist;
else
    cent_lookback = index-c_locs_rel(end);
%     min_cent_dist = c_pks(end);
end

h = figure(1);
imshow(frame); hold on;
scatter(curr_cent(1),curr_cent(2),'g+');
scatter(cent(cent_lookback,1),cent(cent_lookback,2),'ro');
pause();
close(h);


end

