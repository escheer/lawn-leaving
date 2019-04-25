function [neighbors_sq_idx, neighbors_lin_idx] = get_pixel_neighbors( bin_img, seed_sq_idx )
%get_pixel_neighbors.m Takes in an image and a seed pixel (expressed in row,col format).
% Returns the pixel indices of all neighboring pixels in square format
% (row,col) and linear indexing. Excludes any neighbors that go outside the
% boundary of the image.

% seed_lin_idx = sub2ind(size(bin_img), seed_sq_idx(2), seed_sq_idx(1)); 

pad_img = padarray(bin_img, [1 1], 1); %put ones around the border of the image to prevent out-of-bounds errors accessing neighbors
M = size(pad_img,1);
%what are the linear indices of these border pixels (make sure not to
%include them in final tally)?
justborder_idx = find(padarray(zeros(size(bin_img)),[1 1],1)); 
%translate the index of seed pixel to the index in the padded image
pad_seed_lin_idx = sub2ind(size(pad_img), seed_sq_idx(:,2)+1, seed_sq_idx(:,1)+1); %have to flip which pixel goes first between matrix and image coordinates and add 1 to account for padding
%Neighbor Offsets
East =          M;
Southeast =     M + 1;
South =         1;
Southwest =    -M + 1;
West =         -M;
Northwest =    -M - 1;
North =        -1;
Northeast =     M - 1;
%
neighbor_offsets = [East, Southeast, South, Southwest, West, Northwest, North, Northeast];

pad_neighbors = bsxfun(@plus, pad_seed_lin_idx, neighbor_offsets)'; %add the neighbor_offsets to the seed pixel to get neighbor indices
pad_neighbors_lin_idx = setdiff(pad_neighbors,justborder_idx); %remove those pixels that were out of bounds

[col,row] = ind2sub(size(pad_img),pad_neighbors_lin_idx);

neighbors_sq_idx = [row-1 col-1]; %remove padding
neighbors_lin_idx = sub2ind(size(bin_img),neighbors_sq_idx(:,2),neighbors_sq_idx(:,1)); %translate to linear indexing using the original image dimensions
end

