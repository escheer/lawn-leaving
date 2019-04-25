function [orig_background, clean_background, level] = getbackground6_prctile(v, level, p, num_gaussian)
%GETBACKGROUND.M Summary of this function goes here
%   Detailed explanation goes here
number_of_frames = v.NumberOfFrames;
bg_frames = floor(linspace(1,number_of_frames,250)); %was 90
nf = length(bg_frames);

% Fill the background matrix with zeros
stack_back = zeros(v.Height, v.Width, nf);
% Stack all the frames together
for i = 1:nf
    frame = imcomplement(im2double(rgb2gray(read(v,bg_frames(i))))); % read in frame ITS DOUBLE
    stack_back(:,:,i) = frame;
end

% Percentile of the stack to obtain the background
orig_background  = prctile(stack_back,p,3);

clean_background = imgaussfilt(orig_background,num_gaussian);

end

