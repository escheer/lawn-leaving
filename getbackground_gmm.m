function [detector, orig_background] = getbackground_gmm( v, region_rounded, outer_boundary_mask)
%GETBACKGROUND_GMM.M This function trains a gaussian mixture model
%foreground detector, which can be used to track a moving animal in the
%subsequent tracking.

number_of_frames = v.NumberOfFrames;
bg_frames = floor(linspace(1,number_of_frames,floor(number_of_frames/20))); %use these frames to train the foreground detector
nf = length(bg_frames);

learningrate = 0.0005;
detector = vision.ForegroundDetector(...
    'AdaptLearningRate',true,...
    'LearningRate',learningrate,...
    'NumTrainingFrames',nf,...
    'MinimumBackgroundRatio',0.1,...
    'NumGaussians',5,...
    'InitialVariance','Auto');

stack_back = zeros(region_rounded(4)+1, region_rounded(3)+1, nf);
% videoPlayer1 = vision.VideoPlayer();
for i = 1:nf %train the foreground detector model
    fr = bg_frames(i);
    frame = imcomplement(im2double(rgb2gray(read(v,fr)))); % read in frame ITS DOUBLE
    masked = frame.*outer_boundary_mask;
    cropped = imcrop(masked,region_rounded);
    stack_back(:,:,i) = cropped;
    fgMask = step(detector, cropped);
%     step(videoPlayer1, fgMask);
end

orig_background = (mean(stack_back,3));

% videoPlayer2 = vision.VideoPlayer();
% for i = 1:1000
% %         fr = bg_frames(i);
%     fr = i;
%     frame = imcomplement(im2double(rgb2gray(read(v,fr)))); % read in frame ITS DOUBLE
%     masked = frame.*outer_boundary_mask;
%     cropped = imcrop(masked,region_rounded);
%     fgMask = step(detector, cropped);
%     step(videoPlayer2, fgMask);
% end


end

