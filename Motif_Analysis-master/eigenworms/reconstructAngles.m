function reconstructedAngle = reconstructAngles(timeSeries, eigenWorms, ...
    numDimensions)

% RECONSTRUCTANGLES  Use eigenworms to reconstruct angle data.  Save the
% output in the same location as the projected amplitudes file.
%
% Input
%   timeSeries         - The projected amplitudes time series to be
%                        analysed
%   eigenWorms         - the basis eigenWorms that were used in the
%                        decomposition
%   numDimensions      - The number of eigen worms to use for the
%                        reconstruction
%
% Output
%   reconstructedAngle - The array of reconstructed angles for each frame
%                        (dimensions are angleNumber * frameNumber)
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

reconstructedAngle = zeros(size(eigenWorms, 2), size(timeSeries, 2));
%reconstruct worm skeleton angles (add up eigenworms times amplitudes)
for j = 1:size(timeSeries, 2)
    for k = 1:numDimensions
        % multiply the amplitude by each eigenworm and add up the result
        a = timeSeries(k, j) .* eigenWorms(k, :);
        reconstructedAngle(:, j) = reconstructedAngle(:, j) + a';
    end
end