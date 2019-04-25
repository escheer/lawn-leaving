function [spline, spline_wc, curvature] = getWormSpline( wormcrop, x_offset, y_offset )
%GETWORMSPLINE.m this function is only a wrapper for Navin's package for
%body_contour_from_image

body_contour_struct = body_contour_from_image(wormcrop);
spline = [body_contour_struct.x' body_contour_struct.y']; %52 points along the cropped worm image
%trim spline until it is all contained within worm body
if ~isempty(spline) && sum(sum(spline == 0)) ~= numel(spline) %apparently this can happen, idk why...navin
    spline_wc = int32(spline);
    while wormcrop(spline_wc(1,2),spline_wc(1,1))==0
        spline_wc = spline_wc(2:end,:);
        disp('TRIM! from end 1');
    end
    while wormcrop(spline_wc(end,2),spline_wc(end,1))==0
        spline_wc = spline_wc(1:end-1,:);
        disp('TRIM! from end 2');
    end
    spline_wc = double(spline_wc);
    spline = [spline_wc(:,1)+double(x_offset) spline_wc(:,2)+double(y_offset)];
    
    %scale the smoothing #frames by the size of the spline
    framelen = ceil(size(spline_wc,1)/3.4667);
    if mod(framelen,2)==0 %make sure it's odd for the savitzy-golay filter
        framelen = framelen-1;
    end
    smooth_spline = sgolayfilt(double(spline_wc),3,framelen);
    Vertices=smooth_spline;
    Lines=[(1:size(Vertices,1))' (2:size(Vertices,1)+1)']; Lines(end,2)=1;
    curvature=LineCurvature2D(Vertices,Lines); %gets the radius of the tangent circle at every point
else
    spline = NaN;
    spline_wc = NaN;
    curvature = NaN;
end
end

