%check for bumps
clear; clc; close all;
[~, pathname, ~] = uigetfile({'*'}); %this is nice because you can see which ones are in progress
cd(pathname);
movielist = dir('*.avi');
movienames = {movielist.name}';

%     Get next file and delete it from the list of movies to be processed
[ movienames, videoname ] = getnextvideo( movienames );
v = VideoReader([pathname videoname]);
number_of_frames = v.NumberOfFrames;

v2 = VideoReader([pathname videoname]); %have to re-instantiate this object so we can use the faster readFrame method.
v2.CurrentTime = 0; % start at the beginning movie THIS DOESNT EXIST IN 2014a

prev_frame = [];
frame = [];
totalframe = 1;
videoframe = 1;
vidindex = 1;
batch_reg_info(1).videoname = videoname;
offsets = zeros(number_of_frames,6); %column 5 is videoframe, column 6 is totalframe
offsets(1,5)=1; offsets(1,6)=1;%initialize with first frame
tic;
while hasFrame(v2) || ~isempty(movienames) %MAIN LOOP
    if ~hasFrame(v2) %change the video here
        batch_reg_info(vidindex).offsets = offsets; %save the info for the last video
        disp('CHANGE VIDEO!');
        videoframe = 1;
        vidindex = vidindex+1;
        [ movienames, videoname ] = getnextvideo( movienames );
        clear v; clear v2;
        v = VideoReader([pathname videoname]);
        v2 = VideoReader([pathname videoname]);%for tracking
        v2.CurrentTime = 0;
        number_of_frames = v.NumberOfFrames;
        offsets = zeros(number_of_frames,6);
        batch_reg_info(vidindex).videoname = videoname;
    end
    disp(['curr frame is ' num2str(totalframe)]);
    prev_frame = frame;
    frame = rgb2gray(readFrame(v2));
    if totalframe == 1
        x = figure(1); imshow(frame);
        h = imrect(); wait(h);
        crop_pos = getPosition(h);
        frame = imcrop(frame,crop_pos);
        close(x);
    end
    if totalframe>1
        frame = imcrop(frame,crop_pos);
        usfac = 100;
        [output, ~] = dftregistration(fft2(prev_frame),fft2(frame),usfac);
        offsets(videoframe,1:4)=output;
        offsets(videoframe,5)=videoframe;
        offsets(videoframe,6)=totalframe;
    end
    totalframe = totalframe+1;
    videoframe = videoframe+1;
end
batch_reg_info(vidindex).offsets = offsets;
toc;