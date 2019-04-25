function plot_tracks_oneatatime_060718( bg_struct, allTracks )
%PLOT_TRACKS_ONEATATIME.M Simple function that plots each track's head
%position one track at a time.
h1 = figure(1); hold on;
plot(bg_struct(1).ev_ho_crp_rel(:,1),bg_struct(1).ev_ho_crp_rel(:,2));
axis equal;
% h2 = figure(2); hold on;
% plot(bg_struct(1).ev_ho_crp_rel(:,1),bg_struct(1).ev_ho_crp_rel(:,2));
% axis equal;
for i = 1:length(allTracks)
    disp(i);
    figure(h1);
    plot(allTracks(i).centroid(:,1),allTracks(i).centroid(:,2));
%     figure(h2);
%     plot(allTracks(i).head(:,1),allTracks(i).head(:,2));
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
% close(h2);

end

