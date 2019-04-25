function [spline, curvature, posture_angle] = getWormSpline2( worm, x_offset, y_offset )
%GETWORMSPLINE2.m this function calculates the smoothed spline and
%curvature based on the midline pixels identified during segmentation. then
%comes up with curvature

roughspline = worm.skeleton.pixels;
rough_x = roughspline(:,2)+x_offset;
rough_y = roughspline(:,1)+y_offset;
xx = linspace(1,length(rough_x),49);
p = 0.6;
pp_x = csaps(1:length(rough_x),rough_x,p,xx,[5 ones(1,length(rough_x)-2) 5]);
pp_y = csaps(1:length(rough_y),rough_y,p,xx,[5 ones(1,length(rough_y)-2) 5]);

spline = [pp_x' pp_y'];
[posture_angle,~] = angleFromXY(pp_x, pp_y); %AEX code

Vertices = spline;
Lines=[(1:size(Vertices,1))' (2:size(Vertices,1)+1)']; Lines(end,2)=1;
curvature=LineCurvature2D(Vertices,Lines);%gets the radius of the tangent circle at every point
end

