function [minDist bestMatchStart] = motifDist(timeSeries, motif)

% MOTIFDIST Given a time series and a shapelet, find the minimum euclidean
% distance between the time series and the shapelet and the location in the
% time series where the minimum distance occurs.
% 
% Input
%   timeSeries     - the time series to be compared to the motif (must be
%                    as long or longer than motif)
%   motif          - the short time series motif to be compared to the time
%                    series
%   
% Output
%   minDist        - the minimum distance between the motif and the time
%                    series
%   bestMatchStart - the index where the best match occurred
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


% check that timeSeries is as long or longer than the motif
if size(timeSeries, 1) < size(motif, 1)
    error('The time series must be longer than the motif.')
end

motifLength = size(motif, 2);
minDistSquare = Inf;

for i = 1:(size(timeSeries, 2) - motifLength + 1)
    %get the current subsequence
    subsequence = timeSeries(:, i:(i+motifLength-1));
        
    distSquare = sum(sum((subsequence - motif).^2));
    
    if distSquare < minDistSquare
        % if the current distance is less than the previous minimum, update 
        % minDist and bestMatchStart
        minDistSquare = distSquare;
        bestMatchStart = i;
    end
end

minDist = sqrt(minDistSquare);