
%% Update Unassigned Tracks
% Mark each unassigned track as invisible, and increase its age by 1.

% Elias Scheer
% 1-12-17
% borrowed from MultiObjectTracking example

function [tracks, tracksThatLeft] = updateUnassignedTracks_wholevideo(tracks,tracksThatLeft, unassignedTracks,curr_frame)
for i = 1:length(unassignedTracks)
    ind = unassignedTracks(i);
    tracks(ind).age = tracks(ind).age + 1;
    tracks(ind).consecutiveInvisibleCount = tracks(ind).consecutiveInvisibleCount + 1;
    tracks(ind).framesActive = [tracks(ind).framesActive curr_frame]; %by convention the last frame is the frame in which the track was lost.
end
[tracks, tracksThatLeft] = identifyTracksStillIn_wholevideo(tracks, tracksThatLeft);
end

