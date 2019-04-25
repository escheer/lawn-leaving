function [ ev_ho_x, ev_ho_y ] = geteventhorizon(background)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
BW = edge(background,'canny');
BW2 = bwareaopen(BW,50); %get rid of schmutz
[edge_y, edge_x] = find(BW2);

imshow(BW2); hold on;
%select the region in which to find the boundary
H = imfreehand;
Position = wait(H);
pos_x = Position(:,1); pos_y = Position(:,2);

%find those edge pixels which fall inside the polygon specified by
%imfreehand
in = inpolygon(edge_x,edge_y,pos_x,pos_y);
scatter(edge_x(in),edge_y(in),10,'b'); %these points are inside the drawn polygon;
inside_x = edge_x(in);
inside_y = edge_y(in);

%calculate the convex hull of these points to get the lawn boundary.
K = convhull(inside_x,inside_y);
scatter(inside_x(K), inside_y(K),25,'r+');

ev_ho_x = inside_x(K); ev_ho_y = inside_y(K);
close all;
end
