function mat=computeRotationMatrix(theta)
%mat=computeRotationMatrix(theta)
%
% compute 2D rotation matrix for an angle
%
%(c) Saul Kato, saul@kato.com
%Created 2014-3-31

    %this is a counterclockwise rotation
    mat=[cos(theta) -sin(theta); sin(theta) cos(theta) ]';

end