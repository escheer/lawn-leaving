function [e1_call, e2_call, loop1_ctr, loop2_ctr] = get_worm_head(endpt_1,endpt_2,wormcrop_bin,wormcrop_orig)
%get_worm_head.m This function takes in a masked worm and a spline and
%identifies which side is the head by comparing which endpoint of the
%spline is brighter.
worm_area = sum(sum(wormcrop_bin));
threshold = worm_area/9;

% sum the pixel values starting from either end of the worm
% first, find a group of pixels inside the thresholded worm on either end
% of the spline

inv_worm = imcomplement(wormcrop_bin); %worm is in black (0), outside is white(1)

[~,IDX] = bwdist(wormcrop_bin);  %find nearest non-zero neighbor of the endpoints
[i,j] = ind2sub(size(wormcrop_bin),IDX(endpt_1(2),endpt_1(1)));
endpt_1 = [j i];
[i,j] = ind2sub(size(wormcrop_bin),IDX(endpt_2(2),endpt_2(1)));
endpt_2 = [j i];

%find active pixels for endpt_1
active_pixels_sq = endpt_1;
active_pixels_lin = sub2ind(size(wormcrop_bin),endpt_1(2),endpt_1(1));
grabbed_pixels = [];
loop1_ctr = 0;
while length(grabbed_pixels)<=threshold %search for pixels until you exceed an area threshold
    if loop1_ctr>50
        disp('debug this long loop!');
    end
    % Set the active pixels to 1.
    inv_worm(active_pixels_lin) = 1;
    % The new active pixels list is the set of neighbors of the
    % current list.
    [active_pixels_sq, active_pixels_lin] = get_pixel_neighbors(inv_worm, active_pixels_sq);
    active_pixels_lin = active_pixels_lin(:);
    % Remove from the active_pixels list pixels that are already
    % 1.  This is how we keep the flood fill operation "inside"
    % the white boundaries.
    active_pixels_lin(inv_worm(active_pixels_lin)) = [];
    % Remove duplicates from the list.
    active_pixels_lin = unique(active_pixels_lin);
    %makes sure the two sets match based on the linear indexing set
    [tmp_col, tmp_row] = ind2sub(size(inv_worm),active_pixels_lin);
    active_pixels_sq = [tmp_row tmp_col];
    grabbed_pixels = [grabbed_pixels; active_pixels_lin];
    loop1_ctr = loop1_ctr+1;
end
endpt_1_avg = mean(wormcrop_orig(grabbed_pixels));%get average value over active pixels near endpt 1
% top_pixels = sort(wormcrop_orig(grabbed_pixels),'descend');
% endpt_1_avg = median(top_pixels(1:floor(length(top_pixels)/2)));

inv_worm = imcomplement(wormcrop_bin);
%find active pixels for endpt_END
active_pixels_sq = endpt_2;
active_pixels_lin = sub2ind(size(wormcrop_bin),endpt_2(2),endpt_2(1));
grabbed_pixels = [];
loop2_ctr = 0;
while length(grabbed_pixels)<=threshold %search for pixels until you exceed an area threshold
    if loop2_ctr>50
        disp('debug this long loop!');
    end
    % Set the active pixels to 1.
    inv_worm(active_pixels_lin) = 1;
    % The new active pixels list is the set of neighbors of the
    % current list.
    [active_pixels_sq, active_pixels_lin] = get_pixel_neighbors(inv_worm, active_pixels_sq);
    active_pixels_lin = active_pixels_lin(:);
    % Remove from the active_pixels list pixels that are already
    % 1.  This is how we keep the flood fill operation "inside"
    % the white boundaries.
    active_pixels_lin(inv_worm(active_pixels_lin)) = [];
    % Remove duplicates from the list.
    active_pixels_lin = unique(active_pixels_lin);
    %makes sure the two sets match based on the linear indexing set
    [tmp_col, tmp_row] = ind2sub(size(inv_worm),active_pixels_lin);
    active_pixels_sq = [tmp_row tmp_col];
    grabbed_pixels = [grabbed_pixels; active_pixels_lin];
    loop2_ctr = loop2_ctr+1;
end
endpt_2_avg = mean(wormcrop_orig(grabbed_pixels));%get average value over active pixels near endpt 1
% top_pixels = sort(wormcrop_orig(grabbed_pixels),'descend');
% endpt_END_avg = median(top_pixels(1:floor(length(top_pixels)/2)));

% disp(['end 1: ' num2str(endpt_1_avg) '; endpt END: ' num2str(endpt_END_avg)]);
if (endpt_1_avg-endpt_2_avg)>0.005 %if the avg intensity of endpoint 1 is greater, it is the head
    e1_call = 1; %1 means head, 2 means tail
    e2_call = 2;
elseif (endpt_2_avg-endpt_1_avg)>0.005 %vice versa
    e1_call = 2;
    e2_call = 1;
else% if equal, return the endpoint of the last head
    e1_call = NaN;
    e2_call = NaN;
end
end

