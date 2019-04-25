function plot_tracks_oneatatime( allTracks )
%PLOT_TRACKS_ONEATATIME.M Simple function that plots each track's head
%position one track at a time.
h1 = figure(1); hold on;
plot(allTracks(end).ev_ho_x,allTracks(end).ev_ho_y);
axis equal;
h2 = figure(2); hold on;
plot(allTracks(end).ev_ho_x,allTracks(end).ev_ho_y);
axis equal;
for i = 1:length(allTracks)
    disp(i);
    figure(h1);
    plot(allTracks(i).centroid(:,1),allTracks(i).centroid(:,2));
    figure(h2);
    plot(allTracks(i).head(:,1),allTracks(i).head(:,2));
    inp = input('press return for next track, 0 to end!');
    if isempty(inp)
        continue;
    elseif inp == 0
        break;
    else
        continue;
    end
end
close(h1);
close(h2);

end

