function [ min_head_lookback, min_head_dist, min_cent_lookback, min_cent_dist ] = getprevdist( track, index, numlookback, mindist )
%GETPREVDIST.M Given a track, get the distance of the current position
%from the previous locations
head = track.head;
cent = track.centroid;
curr_h = head(index,:);
curr_c = cent(index,:);
head_select = head(index-numlookback:index,:);
cent_select = cent(index-numlookback:index,:);
head_dist = diag(pdist2(head_select,repmat(curr_h,length(head_select),1)));
cent_dist = diag(pdist2(cent_select,repmat(curr_c,length(cent_select),1)));
%find peaks
[h_pks,h_locs] = findpeaks(head_dist,'MinPeakHeight',mindist,'MinPeakProminence',5);
h_locs_rel = h_locs+(index-numlookback);
[c_pks,c_locs] = findpeaks(cent_dist,'MinPeakHeight',mindist,'MinPeakProminence',5);
c_locs_rel = c_locs+(index-numlookback);

if isempty(h_pks)
%     [~,head_closest] = nanmin(abs(head_dist-(mindist+5)));
    [~,head_closest] = nanmax(abs(head_dist-(mindist+5))<=1);
    if head_dist(head_closest)<mindist %if the best you can do is less than mindist+5, just take the max distance
        [~, head_closest] = nanmax(head_dist);
    end
    head_closest_rel = head_closest+(index-numlookback);
    min_head_lookback = index-head_closest_rel;
    min_head_dist = head_dist(head_closest);
    h_locs = head_closest;
    h_pks = min_head_dist;
else
    min_head_lookback = index-h_locs_rel(end);
    min_head_dist = h_pks(end);
end

if isempty(c_pks) %either all the peaks are too low or its just smoothly decreasing i.e. when the worm is moving fast
%     [~,cent_closest] = nanmin(abs(cent_dist-(mindist+5)));
    [~,cent_closest] = nanmax(abs(cent_dist-(mindist+5))<=1);
    if cent_dist(cent_closest)<mindist %if the best you can do is less than mindist+5, just take the max distance
        [~, cent_closest] = nanmax(cent_dist);
    end
    cent_closest_rel = cent_closest+(index-numlookback);
    min_cent_lookback = index-cent_closest_rel;
    min_cent_dist = cent_dist(cent_closest);
    c_locs = cent_closest;
    c_pks = min_cent_dist;
else
    min_cent_lookback = index-c_locs_rel(end);
    min_cent_dist = c_pks(end);
end

% h = figure(); hold on;
% plot(cent_dist); scatter(c_locs,c_pks,'k+');
% plot(head_dist); scatter(h_locs,h_pks,'r+');
% legend('cent dist','cent peaks','head dist','head peaks');

end

