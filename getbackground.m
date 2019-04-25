function [orig_background, clean_background, bgthresh] = getbackground( v, region_rounded, outer_boundary_mask, outer_boundary_line, old_bgthresh, skip_choose )
%GETBACKGROUND.M Summary of this function goes here
%   Detailed explanation goes here
number_of_frames = v.NumberOfFrames;
bg_frames = floor(linspace(1,number_of_frames,90));
nf = length(bg_frames);
% Fill the background matrix with zeros
% stack_back = zeros(region_rounded(4)+1, region_rounded(3)+1, nf);
stack_back = zeros(v.Height, v.Width, nf);
% Stack all the frames together
for i = 1:nf
    frame = imcomplement(im2double(rgb2gray(read(v,bg_frames(i))))); % read in frame ITS DOUBLE
    stack_back(:,:,i) = frame;
end

% Mean of the stack to obtain the background
orig_background = (mean(stack_back,3));

bgcopy = orig_background;
if ~skip_choose
    thresh = multithresh(bgcopy,2);
    [~, bgthresh, ~] = choose_thresh(bgcopy,thresh(2));
    close all;
    abovethresh = im2bw(orig_background,bgthresh);
    % % don't inpaint NaNs from the outer black color
    [abvpix_y, abvpix_x] = find(abovethresh);
    dist = pdist2([abvpix_x abvpix_y],outer_boundary_line);
    [revert_i, ~] = find(dist<100);
    abovethresh(sub2ind(size(abovethresh),abvpix_y(revert_i),abvpix_x(revert_i)))=0;
    % %
    bgcopy(abovethresh) = NaN;
    filledbg = inpaint_nans(bgcopy);
    figure(); subplot(121);
    imshow(orig_background);
    subplot(122);
    imshow(filledbg);
    filled_or_not = input('Which is better background LEFT (1) or RIGHT (2)?');
    if filled_or_not == 2
        clean_background = filledbg;
%         clean_is_edited = 1;
    elseif filled_or_not == 1
        orig_background = orig_background.*outer_boundary_mask; %apply mask
        orig_background = imcrop(orig_background,region_rounded); %crop it
        clean_background = orig_background;
        bgthresh = NaN;
%         clean_is_edited = 0;
        close all;
        return;
    else
        error(' you must choose one or 2!');
    end
    close all;
else
    bgthresh = old_bgthresh;
    if ~isnan(bgthresh) %this means that on the first go-round, you chose to keep the backround without cleaning it!
        abovethresh = im2bw(orig_background,bgthresh);
        % % don't inpaint NaNs from the outer black color
        [abvpix_y, abvpix_x] = find(abovethresh);
        dist = pdist2([abvpix_x abvpix_y],outer_boundary_line);
        [revert_i, ~] = find(dist<100);
        abovethresh(sub2ind(size(abovethresh),abvpix_y(revert_i),abvpix_x(revert_i)))=0;
        % %
        bgcopy(abovethresh) = NaN;
        filledbg = inpaint_nans(bgcopy);
        clean_background = filledbg;
    else
        clean_background = orig_background;
    end
end

clean_background = imgaussfilt(clean_background,10);%blur out any fine features of the background with a radius much larger than the worm
clean_background = clean_background.*outer_boundary_mask; %apply mask
clean_background = imcrop(clean_background,region_rounded); %crop it

orig_background = orig_background.*outer_boundary_mask; %apply mask
orig_background = imcrop(orig_background,region_rounded); %crop it


end

