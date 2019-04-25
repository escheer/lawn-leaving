function [kappa, xx, yy] = relative_curvature_angles_from_coordinates(input_x, input_y, input_num_contour_points)
% kappa = relative_curvature_angles_from_coordinates(input_x, input_y, input_num_contour_points)

if(length(input_x)<=5)
    kappa = zeros(1,input_num_contour_points);
    xx = kappa;
    yy = kappa;
    return;
end

num_contour_points = input_num_contour_points+2;

if ismac
    x = smooth(input_x); %this works with the curve fitting 
    y = smooth(input_y);
else
    x = movmean(input_x,5);%changed 08/07/17 by Elias Scheer to avoid using the Curve Fitting Toolbox
    y = movmean(input_y,5);    
end

% spline to add more points to contour
t = 1:length(x);
ts = 1:(length(x)/(num_contour_points-1)):(length(x)+1);
if(ts(end)<length(x))
    ts = [ts length(x)];
end
xx = (spline(t,x,ts));
yy = (spline(t,y,ts));

kappa = zeros(1,length(xx))+NaN;
for(k=2:length(xx)-1)
    kappa(k) = ((180/pi)*atan2(-(yy(k-1)-yy(k)),(xx(k-1)-xx(k))));
end
dkappa = abs([NaN diff(kappa)]);
jump_idx = find((dkappa)>=90);
original_jump_len = length(jump_idx);
ctr=0;
while(~isempty(jump_idx))
    idx = jump_idx(1);
    if(idx>1)
        if(kappa(idx-1)>0)
            kappa(idx) = kappa(idx) + 360;
        else
            kappa(idx) = -360 + kappa(idx);
        end
    end
    dkappa = abs([NaN diff(kappa)]);
    jump_idx = find((dkappa)>=90);
    ctr = ctr + 1;
    if(ctr>1000*original_jump_len)
        break;
    end
end
kappa(find(kappa>360)) = 360;
kappa(find(kappa<-360)) = -360;
kappa = kappa - nanmedian(kappa);
kappa(find(kappa>90)) = 90;
kappa(find(kappa<-90)) = -90;

kappa = matrix_replace(kappa,'==',NaN,0); 
kappa(1)=[]; 
kappa(end)=[]; 

return;
end


%         else
%             % curvature
%             % derivatives at each point
%             xp = smooth([0; diff(xx)']);
%             yp = smooth([0; diff(yy)']);
%             xpp = smooth([0; diff(xp)]);
%             ypp = smooth([0; diff(yp)]);
%             kappa = (xp.*ypp)-(yp.*xpp) ./ (xp .^2 + yp .^2) .^ (3/2);
%             kappa(find(kappa>0.5)) = 0.5;
%             kappa(find(kappa<-0.5)) = -0.5;
%             kappa(find(isnan(kappa))) = 0;
%             kappa = kappa/Track.PixelSize; % kappa is kappa = 1/pixel; convert to 1/mm
%         end
