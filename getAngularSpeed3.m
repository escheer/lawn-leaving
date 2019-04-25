function [AngSpeed, AngSpeed_rad] = getAngularSpeed3(centroid,oversample_factor)
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

if oversample_factor == 1
    Vertices = centroid;
else
    xx = linspace(1,length(cent_x),length(cent_x)*oversample_factor);
    p = 1; %perfect adherence to data
    pp_x = csaps(1:length(cent_x),cent_x,p,xx,[5 ones(1,length(cent_x)-2) 5]); %cubic spline interpolation
    pp_y = csaps(1:length(cent_y),cent_y,p,xx,[5 ones(1,length(cent_x)-2) 5]);
    Vertices = [pp_x' pp_y'];
end

%approach 1
Lines= [(1:(size(Vertices,1)-1))' (2:size(Vertices,1))'];
curvs = LineCurvature2D(Vertices,Lines); %get curvature at each point k = d_theta/d_s
d_s = sqrt( sum( abs( diff( Vertices ) ).^2, 2 ) ); %d_s = arc length from one point to the next
d_theta = curvs(2:end,:).*d_s; %arc length * k = d_theta (ang speed)
d_theta(d_theta>pi)=pi; %maximum angular speed is pi or -pi
d_theta(d_theta<-pi)=-pi;

if oversample_factor ==1
    as_down = [NaN;d_theta];
else
    [pos_pks,pos_locs] = findpeaks(d_theta,'MinPeakHeight',0.05); %make sure you don't lose peak height after resizing
    [neg_pks,neg_locs] = findpeaks(-1*d_theta,'MinPeakHeight',0.05);
    as_down = imresize(d_theta,[size(centroid,1) 1],'Bilinear'); %resize vector by averaging

    pos_locs_down = round(pos_locs./oversample_factor);
    neg_locs_down = round(neg_locs./oversample_factor);
    %make sure there are no zero indices
    pos_locs_down(pos_locs_down==0)=1;
    neg_locs_down(neg_locs_down==0)=1;
    
    [pos_locs_down, pl_idx] = unique(pos_locs_down);
    [neg_locs_down, nl_idx] = unique(neg_locs_down);
    pos_pks_down = pos_pks(pl_idx);
    neg_pks_down = neg_pks(nl_idx);
    as_down(pos_locs_down) = pos_pks_down;
    as_down(neg_locs_down) = -1*neg_pks_down;
end

AngSpeed_rad = as_down;
AngSpeed = -180/pi.*AngSpeed_rad;

end

 