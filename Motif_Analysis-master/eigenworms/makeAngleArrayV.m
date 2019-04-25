function [angleArray, meanAngles] = makeAngleArrayV(x, y)

%MAKEANGLEARRAY Get tangent angles for each frame of normBlocks and rotate
%               to have zero mean angle
%
%   [ANGLEARRAY, MEANANGLES] = MAKEANGLEARRAY(X, Y)
%
%   Input:
%       x - the x coordinates of the worm skeleton (equivalent to
%           dataBlock{4}(:,1,:)
%       y - the y coordinates of the worm skeleton (equivalent to
%           dataBlock{4}(:,2,:)
%
%   Output:
%       angleArray - a numFrames by numSkelPoints - 1 array of tangent
%                    angles rotated to have mean angle of zero.
%       meanAngles - the average angle that was subtracted for each frame
%                    of the video.
%
% Copyright Medical Research Council 2013
% Andr� Brown, andre.brown@csc.mrc.ac.uk, aexbrown@gmail.com
%
%
% The MIT License
%
% Copyright (c)  Medical Research Council 2013
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.


lengthX = size(x,2);

% calculate the x and y differences
dX = diff(x,1,2);
dY = diff(y,1,2);

% calculate tangent angles.  atan2 uses angles from -pi to pi instead...
% of atan which uses the range -pi/2 to pi/2.
angles = atan2(dY, dX);

% need to deal with cases where angle changes discontinuously from -pi
% to pi and pi to -pi.  
angles = unwrap(angles,[],2);

% rotate skeleton angles so that mean orientation is zero
meanAngles = mean(angles,2);
angles = angles - meanAngles(:,ones(1,lengthX-1));

angleArray = angles;
