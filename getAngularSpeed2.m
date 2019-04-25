function [AngSpeed, AngSpeed_rad] = getAngularSpeed2(centroid, interval)
%GETANGULARSPEED2.m extracts the angular speed in deg/sec from the centroid
%path (there should be n-2 measurements for this for a centroids vector of
%length n.
AngSpeed = NaN(size(centroid,1),1);
AngSpeed_rad = NaN(size(AngSpeed));

centroid_os1 = centroid(1+interval:end,:);
centroid_os2 = centroid(1+2*interval:end,:);

for idx = 1:size(centroid_os2,1)

    % a,b,c are coordinates of 3 time points
    a = centroid(idx,:);
    b = centroid_os1(idx,:);
    c = centroid_os2(idx,:);
    
%     DirVector1=[b-a 0];
%     DirVector2=[c-b 0];

    DirVector3 = b-a;
    DirVector4 = c-b;
    AngSpeed(idx+interval) = -180/pi*atan2( DirVector4(1)*DirVector3(2)-DirVector4(2)*DirVector3(1) , DirVector4(1)*DirVector3(1)+DirVector4(2)*DirVector3(2));
    AngSpeed_rad(idx+interval) = atan2( DirVector4(1)*DirVector3(2)-DirVector4(2)*DirVector3(1) , DirVector4(1)*DirVector3(1)+DirVector4(2)*DirVector3(2));
    
    % this method doesn't give a sign to the angle.
    %     AngSpeed2(idx+3) = 180/pi*(atan2(norm(cross(DirVector1,DirVector2)),dot(DirVector1,DirVector2)));
end

end

 