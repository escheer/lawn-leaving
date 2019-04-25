function featureMatrix = motifFeatureMatrix(fileList, motifDictionary, ...
    calculateFreq, numDimensions, samplingInterval, lengthThreshold, display)

% MOTIFFEATUREMATRIX Find the minimum distance, best match location, and
% (optionally) the frequency of matches between each file in file list and
% each motif in motifDictionary.  Each motif is compared against every
% subsequence of the same length in each time series file in fileList.
% Euclidean distance is used in all cases.  The output is a feature matrix
% with time series files along the rows and distances to motifs along the
% columns.
%
% Input
%   fileList         - List of the full paths of the files to be processed.
%   motifDictionary  - A cell array containing the parent file names, the
%                      actual motifs, and their start locations at the
%                      original time series resolution (i.e. no longer down
%                      sampled)
%   calculateFreq    - Logical.  Should the frequency of motif occurances
%                      be calculated?  This is the number of times a motif
%                      is present in a time series within a threshold
%                      currently given by 4 time the motif length.
%   numDimensions    - The number of eigen worms to consider in motif
%                      finding
%   samplingInterval - Used to downsample data for faster processing. Eg. a
%                      value of 6 means take every sixth point from the
%                      timeseries.
%   lengthThreshold  - Time series shorter than this value will not be
%                      included in the analysis.
%   display          - 0, 1, or 2.  Should the progress be displayed in the
%                      command window? 0 = no, 1 = display progress per
%                      file, 2 = display more detailed progress (file,
%                      current motif from dictionary)
%
% Output
%   featureMatrix    - A matrix of distances (and frequencies of
%                      occurances) between each file in file list and each
%                      motif in motifDictionary
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


% initialize a feature matrix (must be twice as large if frequencies are to
% be calculated)
if calculateFreq
    featureMatrix = NaN(size(fileList, 1), size(motifDictionary, 1) * 2);
else
    featureMatrix = NaN(size(fileList, 1), size(motifDictionary, 1));
end

% loop through fileList
for i = 1:size(fileList, 1)
    % display progress
    if display == 1
        disp(i/size(fileList, 1))
    end
    
    % load the time series to examine
    timeSeries = cell2mat(struct2cell(load( fileList{i} )));
    
    % check that timeSeries is not empty and not all NaNs.  We also don't
    % want any data significantly shorter than 15 minutes
    if length(timeSeries) > lengthThreshold && any(~isnan(timeSeries(:)))
        
        % downsample the amplitude data for comparison with motif
        timeSeriesDownSamp = ...
            timeSeries(1:numDimensions, 1:samplingInterval:end);
        
        % calculate the distance from each of the motifs
        for j = 1:size(motifDictionary, 1)
            % display more detailed progress
            if display == 2
                disp([num2str(i/size(fileList, 1)) ', ' ...
                    num2str(j/size(motifDictionary, 1))])
            end
            % check that motifDictionary entry is not empty
            if ~isempty(motifDictionary{j, 1})
                
                % calculate the distance from one of the motifs.  They
                % should be fairly similar by definition and this saves
                % time compared to calculated both distances and taking the
                % minimum.
                minDist = motifDistance(timeSeriesDownSamp', ...
                    size(timeSeriesDownSamp, 2), numDimensions,...
                    motifDictionary{j, 5}(1:numDimensions, ...
                    1:samplingInterval:end)', ...
                    size(motifDictionary{j, 5}(1:numDimensions, ...
                    1:samplingInterval:end), 2));
                                
                % normalise distance by the motif length
                minDistNorm = minDist/sqrt(size(motifDictionary{j, 3}, 2));
                
                % add the distance to the shapelet feature matrix
                featureMatrix(i, j) = minDistNorm;
                
                % if requested, calculate the frequency of occurrance
                if calculateFreq
                    % define the distance threshold based on the motif length
                    matchThreshold = ...
                        0.01*size(motifDictionary{j, 3}(1:numDimensions,...
                        1:samplingInterval:end), 2);
                    
                    % Because both motifs in the motif dictionary should be
                    % quite similar and we're allowing significant slack in
                    % finding a match, just use the first motif.
                    freq = motifFrequency(timeSeriesDownSamp,...
                        motifDictionary{j, 3}(1:numDimensions,...
                        1:samplingInterval:end), matchThreshold);
                    
                    % add the frequencies after the distances
                    featureMatrix(i, j + size(motifDictionary, 1)) = freq;
                end
            end
        end
    end
end
