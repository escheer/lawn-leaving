function [orig_background, clean_background, filledbg, level] = getbackground4_prctile(v, level, p)
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
% orig_background = mean(stack_back,3);OLDWAY

% clean_background = orig_background;
% abovethresh = zeros(size(clean_background));
% dwellingworms = zeros(size(abovethresh));

clean_background = imgaussfilt(orig_background,10);

%find any pixels above threshold in the lawn to fill in and then blur the
%image with a gaussian filter to generate the clean background
% cln_bgsub_stack = zeros(v.Height, v.Width, nf);
% for i = 1:nf
%     cln_bgsub_stack(:,:,i) = imcomplement(stack_back(:,:,i)-clean_background);
% end
% avg_clean_bgsub = mean(cln_bgsub_stack,3);
% 
% abovethresh = imcomplement(im2bw(avg_clean_bgsub,level));
% dwellingworms = zeros(size(abovethresh));
% if firstvideo
%     figure();
%     ha = tight_subplot(1,2);
%     axes(ha(1)); imshow(orig_background);
%     axes(ha(2)); imshow(orig_background.*imcomplement(abovethresh));
%     goodthresh = input('BACKGROUND SCHMUTZ WELL-THRESHOLDED? (1/0)');
%     while ~(goodthresh==0 || goodthresh == 1)
%         goodthresh = input('BACKGROUND SCHMUTZ WELL-THRESHOLDED? (1/0)');
%     end
%     close all;
%     if ~goodthresh
%         [~, level, ~] = choose_thresh(imcomplement(avg_clean_bgsub),1-level);
%         level = 1-level;
%     end
%     abovethresh = imcomplement(im2bw(avg_clean_bgsub,level));
%     figure(); imshow(abovethresh);
%     torestrict = input('NEED TO REMOVE DWELLING WORMS FROM BACKGROUND? (1/0)');
%     while ~(torestrict==0 || torestrict == 1)
%         torestrict = input('NEED TO REMOVE DWELLING WORMS FROM BACKGROUND? (1/0)');
%     end
%     close all;
%     if torestrict
%         figure(); imshow(abovethresh);
%         h = imfreehand(); wait(h); vert = getPosition(h);
%         restrictmask = poly2mask(vert(:,1),vert(:,2),size(abovethresh,1),size(abovethresh,2));
%         dwellingworms = abovethresh.*restrictmask;
%         close all;
%     end
% end
% se = strel('disk',1); %helps to ensure that you don't miss any pixels, and ALL of the dirt gets cleaned up
% abovethresh = imdilate(abovethresh,se);
% abovethresh = imdilate(abovethresh,se);
% 
% [abvpix_y, abvpix_x] = find(abovethresh);
% dist = pdist2([abvpix_x abvpix_y],outer_boundary_line);
% [revert_i, ~] = find(dist<20); %don't flip stuff too close to the outer boundary
% abovethresh(sub2ind(size(abovethresh),abvpix_y(revert_i),abvpix_x(revert_i)))=0;
% 
% % intersect abovethresh with the last abovethresh to ensure that you are
% % only flipping pixels of objects that don't move across all videos (user
% % must be careful not to select any worm dwelling pixels on the first
% % video)
% if ~firstvideo
%     dwellingworms = (abovethresh&~old_abovethresh); %those pixels which appeared newly in this video's abovethresh are usually dwelling worms
%     abovethresh = abovethresh&old_abovethresh;
% end
% %slightly enlarge the thresholded worm to aid in regionfill operation
% se = strel('disk',1);
% dwellingworms = imdilate(dwellingworms,se);
% dwellingworms = imdilate(dwellingworms,se);
% 
% abovethresh = abovethresh&~dwellingworms; %remove dwelling worms from abovethresh

% abovethresh = abovethresh.*outer_boundary_mask; %only flip those pixels inside the region of interest for tracking
% dwellingworms = dwellingworms.*outer_boundary_mask; %same

% filledbg1 = regionfill(orig_background,abovethresh.*outer_boundary_mask);
% filledbg2 = regionfill(filledbg1,dwellingworms.*outer_boundary_mask); %blur out dwelling worms, but don't also blur them out of the actual frames while tracking

filledbg = orig_background;
% 
clean_background = imgaussfilt(filledbg,10);%blur out any fine features of the background with a radius much larger than the worm

% % clean_background = clean_background.*outer_boundary_mask; %apply mask
% % clean_background = imcrop(clean_background,region_rounded); %crop it
% 
% % orig_background = orig_background.*outer_boundary_mask; %apply mask
% % orig_background = imcrop(orig_background,region_rounded); %crop it

end

