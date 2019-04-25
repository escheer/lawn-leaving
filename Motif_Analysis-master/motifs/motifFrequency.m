function matchNumber = motifFrequency(timeSeries, motif, matchThreshold)

% MOTIFFREQUENCY Use brute force search to find how often motif appears in
% timeSeries within the Euclidean distance matchThreshold.  To avoid a
% large number of trivial matches, once a match is called, the search is
% advanced outside of the matching subsequence.
% 
% Input
%   timeSeries     - the time series to be compared to the motif (must be
%                    as long or longer than motif)
%   motif          - the short time series motif to be compared to the time
%                    series
%   matchThreshold - the distance between the motif and a subsequence of
%                    timeSeries must be less than this threshold to call a
%                    match
% 
% Output
%   matchNumber    - the number of times the motif came within
%                    matchThreshold in timeSeries
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




% get the dimensions of the time series
ncol = size(timeSeries, 2);
motifLength = size(motif, 2);

% initialisation
matchNumber = 0;
matchThresholdSquared = matchThreshold^2;

i = 1;
% loop through each subsequence
while i <= ncol - motifLength + 1
    % calculate the distance between the current subsequence and the motif
    distSquare = sum(sum((timeSeries(:, i:(i+motifLength-1)) - motif).^2));
    if distSquare < matchThresholdSquared
        % if the current distance is below one of the thresholds, update
        % its frequency of occurance
        matchNumber = matchNumber + 1;
        
        % to avoid further trivial matches, increment i by
        % motifLength - 1.  The -1 is there because the normal increment
        % follows afterwards anyway.
        i = i + motifLength - 1;
    end
    i = i + 1;
end