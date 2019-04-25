%little script to plot tracks and entering/exiting events on top of a frame
%from the video.

%load video
[filename, pathname, ~] = uigetfile({'*'});
v = VideoReader([pathname filename]);
number_of_frames = v.NumberOfFrames;

%load tracking data
[filename, pathname, ~] = uigetfile({'*.mat'});
load([pathname filename]);

%% plot the tracks
last_frame = im2double(read(v,number_of_frames));
figure();
imshow(imadjust(rgb2gray(last_frame)));
hold on;
plot(ev_ho_x+10,ev_ho_y+10) %get rid of these offsets or save information about how the video was cropped!!!
% for i = 1:length(tracksThatLeft)
%     plot(tracksThatLeft(i).centroid(:,1),tracksThatLeft(i).centroid(:,2));
% end
for i = 1:length(tracks)
    plot(tracks(i).centroid(:,1)+10,tracks(i).centroid(:,2)+10,'Color','k','LineWidth',2);
end

%% plot entering and exiting events
%institute a way to look up which track number based on the track id!
track_id = 7;
track_start = tracks(1).framesActive(1);

enter_frames = enter_events(enter_events(:,2)==track_id,1);
enter_frames = enter_frames - track_start;
exit_frames = exit_events(exit_events(:,2)==track_id,1);
exit_frames = exit_frames - track_start;

scatter(tracks(1).centroid(enter_frames,1)+10,tracks(1).centroid(enter_frames,2)+10,50,'g', 'filled');
scatter(tracks(1).centroid(exit_frames,1)+10,tracks(1).centroid(exit_frames,2)+10,50,'r', 'filled');