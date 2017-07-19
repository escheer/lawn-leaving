% countBlobsInOut.m
% 04/26/2017
% This function takes in the tracks which have been saved over the course
% of the video and outputs a count of the number of blobs found inside and
% outside the event horizon at every frame. Note that the blobs outside the
% event horizon will be much more reliable than those inside since there
% will be many multiblobs inside. The number of blobs outside can be used
% to calculate an aversion ratio if the total number of animals in the
% video is known.
function [blobs_in,blobs_out] = countBlobsInOut(tracks,tracksThatLeft,ev_ho_x, ev_ho_y, last_frame_idx )
x = ev_ho_x';
y = ev_ho_y';

tracksToCheck = [tracks, tracksThatLeft]; %merge tracks and tracksThatLeft

blobs_in = zeros(1,last_frame_idx);
blobs_out = zeros(1,last_frame_idx);

for i = 1:length(tracksToCheck)
    firstframe = tracksToCheck(i).framesActive(1);
    
    bbox = tracksToCheck(i).bbox;
    bbox_verts =   [[bbox(:,1) bbox(:,2)]...          % top left
        [bbox(:,1) bbox(:,2)+bbox(:,4)]...            % bottom left
        [bbox(:,1)+bbox(:,3) bbox(:,2)]...            % top right
        [bbox(:,1)+bbox(:,3) bbox(:,2)+bbox(:,4)]];   % bottom right
    
    topleft_out = ~inpolygon(double(bbox_verts(:,1)),double(bbox_verts(:,2)),x,y);
    bottomleft_out = ~inpolygon(double(bbox_verts(:,3)),double(bbox_verts(:,4)),x,y);
    topright_out = ~inpolygon(double(bbox_verts(:,5)),double(bbox_verts(:,6)),x,y);
    bottomright_out = ~inpolygon(double(bbox_verts(:,7)),double(bbox_verts(:,8)),x,y);
    
    bboxes_out = find(topleft_out&bottomleft_out&topright_out&bottomright_out)+firstframe-1; %frame indices
    bboxes_in = find(~topleft_out&~bottomleft_out&~topright_out&~bottomright_out)+firstframe-1; %frame indices
    
    %increment in count and outcount for every frame on the list
    blobs_in(bboxes_in) = blobs_in(bboxes_in)+1;
    blobs_out(bboxes_out) = blobs_out(bboxes_out)+1;
end

































end