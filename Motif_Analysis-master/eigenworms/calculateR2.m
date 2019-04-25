function rSquared = calculateR2(angleArray, reconstructedAngle)

% RECONSTRUCTIONR2  Given a raw angle array and the eigenworm reconstructed
% version, calculate the R^2 value to summarize the correlation
%
% Input
%   angleArray         - The raw skeleton angle array taken from the
%                        feature file
%   reconstructedAngle - The angle array reconstructed from the eigenworm
%                        decomposition
%
% Output
%   rSquared           - R^2 value for the correlation between the raw 
%                        skeleton angles and the reconstructed angles from 
%                        the eigenworm decomposition.
% 
% 
% Copyright Medical Research Council 2013
% André Brown, andre.brown@csc.mrc.ac.uk, aexbrown@gmail.com
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


% only use the values in angleArray that are not NaN
anglesNoNan = angleArray(~isnan(angleArray));

% get the entries in reconstructedAngle at the corresponding (non-NaN)
% locations
reconAngles = reconstructedAngle(~isnan(angleArray));

% get the correlation coefficient
R = corrcoef(anglesNoNan,reconAngles);

% convert to rSquared
rSquared = R(1,2)^2;