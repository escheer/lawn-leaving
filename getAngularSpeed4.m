function [AngSpeed, AngSpeed_rad] = getAngularSpeed4(centroid,downsample_factor)
%GETANGULARSPEED2.m extracts the angular speed in deg/sec from the centroid
%path. This method smoothly interpolates the centroid path at 5x the
%sampled rate to overcome aliasing problems. Then uses LineCurvature to
%extract theta at every point (distance/tangent circle radius), this can
%then be smoothly differentiated to come up with angular speed.

cent_x = centroid(:,1);
cent_y = centroid(:,2);

%fill in any missing data
cent_x = fillmissing(cent_x,'linear');
cent_y = fillmissing(cent_y,'linear');

%downsample centroid before curvature calculation.
cent_x_dwn = downsample(cent_x,downsample_factor);
cent_y_dwn = downsample(cent_y,downsample_factor);

Vertices = [cent_x_dwn cent_y_dwn];
Lines= [(1:(size(Vertices,1)-1))' (2:size(Vertices,1))'];
curvs = LineCurvature2D(Vertices,Lines); %get curvature at each point k = d_theta/d_s
d_s = sqrt( sum( abs( diff( Vertices ) ).^2, 2 ) ); %d_s = arc length from one point to the next
d_theta = curvs(2:end,:).*d_s; %arc length * k = d_theta (ang speed)
d_theta(d_theta>pi)=pi; %maximum angular speed is pi or -pi
d_theta(d_theta<-pi)=-pi;

%expand d_theta back to size of original centroid vector
[pos_pks,pos_locs] = findpeaks(d_theta,'MinPeakHeight',0.05); %make sure you don't lose peak height after resizing
[neg_pks,neg_locs] = findpeaks(-1*d_theta,'MinPeakHeight',0.05);
as_resize = imresize(d_theta,[size(centroid,1) 1],'nearest'); %resize vector by nearest neighbor


AngSpeed_rad = as_resize;
AngSpeed = -180/pi.*AngSpeed_rad;

end

 