function [tracksStillGood, tracksThatLeft] = identifyTracksStillIn_wholevideo(tracks, tracksThatLeft)
%
%identifyTracksStillIn.m


still_visible = NaN(length(tracks),1);

% figure(3); clf; axis equal; hold on;
for i = 1:length(tracks)
    still_visible(i) = ~(tracks(i).consecutiveInvisibleCount > 2); %was 2 before
end

stillGoodInds = find(still_visible);
tracksStillGood = tracks(stillGoodInds); %these should be tracked


newTracksThatLeft = tracks(setdiff(1:length(tracks),stillGoodInds));
% for i = 1:length(newTracksThatLeft)
%     disp(['a track left!: ' num2str(newTracksThatLeft(i).id)]);
%     if newTracksThatLeft(i).id == 72
%         disp('DEBUG');
%     end
% end
WorthSavingTracksThatLeft_idx = [newTracksThatLeft(:).age]'>50; %only save tracks that left which were at least 50 frames long
% if sum(WorthSavingTracksThatLeft_idx)>=1
%     disp('it was worth saving!');
% end
% WorthSavingTracksThatLeft_idx = ones(length(newTracksThatLeft),1);


newTracksThatLeft = newTracksThatLeft(WorthSavingTracksThatLeft_idx);
tracksThatLeft = [tracksThatLeft, newTracksThatLeft];
end

