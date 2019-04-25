function motifIndex = findMotifsFileList(fileList, numDimensions, ...
    samplingInterval, motifLengths, method, display)

% FINDMOTIFSFILELIST  For each class of worm, get a random sample of
% filesPerClass timeseries files.  For each of these files, extract each
% motif for the lengths defined in motifLengths.  Motifs in this case are
% defined as the most repetitive subsequence within the (multi-dimensional)
% timeseries calculated without normalisation.
%
% Input
%   fileList         - The full path of files to analyse
%   numDimensions    - The number of eigen worms to consider in motif
%                      finding
%   samplingInterval - Used to downsample data for faster processing. Eg. a
%                      value of 6 means take every sixth point from the
%                      timeseries.
%   motifLengths     - A vector of motif lengths to search for.  Eg. motif
%                      lengths [5 10 15] means find motifs with lengths 5,
%                      10, and 15 in each timeseries.
%   method           - String indicating the motif search method to be
%                      used. 'brute force', 'MK', or 'adaptive'.  Both
%                      methods should give identical results but may have
%                      different speed.  If the motif is a significant
%                      fraction of the total length, brute force may be
%                      faster. 'adaptive' attempts to use the faster of the
%                      two based on the motif and time series length.
%   display          - Logical.  Should the progress be displayed in the
%                      command window?
%
% Output
%   motifIndex       - A cell containing the start indices of the motifs,
%                      their lengths, their distance from each other, and
%                      the sampling interval used when finding them. The
%                      first entry of each row contains the file name of
%                      the motif source file.
%                      *Note*: all values will be in terms of the
%                      downsampled data so they must be converted to find
%                      the full length motifs in the original time series.
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


% initialise motifIndex
motifIndex = cell(length(fileList)*length(motifLengths), 5);

% count for going through motifIndex
count = 1;
errorCount = 0;

% go through file list
for i = 1:length(fileList)
    try
        % load the time series
        timeSeries = cell2mat(struct2cell(load( fileList{i} )));
        
        % downsample the timeseries and invert for use in motif finding
        % algorithms
        timeSeries = timeSeries(1:numDimensions, 1:samplingInterval:end)';
        
        %for each file find several motifs with different lengths
        for j = 1:length(motifLengths)
            try
                % check if the progress should be displayed
                if display
                    disp((count + errorCount)...
                        /(length(fileList)*length(motifLengths)))
                end
                
                % find the multidimensional motif without normalisation
                switch method
                    case 'brute force'
                        [motifStart1, motifStart2, minDist] = ...
                            keogh_multi_bf(timeSeries, ...
                            size(timeSeries, 1), ...
                            size(timeSeries, 2), motifLengths(j), ...
                            motifLengths(j), 0);
                        
                    case 'MK'
                        [motifStart1, motifStart2, minDist] = ...
                            keogh_multi_mk(timeSeries, ...
                            ones(size(timeSeries,1), 1), ...
                            size(timeSeries, 1), ...
                            size(timeSeries, 2), motifLengths(j), 5, ...
                            motifLengths(j), 0);
                        
                    case 'adaptive'
                        % if the motif is more than 1.5% the length of the
                        % time series, use the brute force method,
                        % otherwise, use MK method.
                        if motifLengths(j)/size(timeSeries, 2) > 0.015
                            [motifStart1, motifStart2, minDist] = ...
                                keogh_multi_bf(timeSeries, ...
                                size(timeSeries, 1), ...
                                size(timeSeries, 2), motifLengths(j), ...
                                motifLengths(j), 0);
                        else
                            [motifStart1, motifStart2, minDist] = ...
                                keogh_multi_mk(timeSeries, ...
                                ones(size(timeSeries,1), 1), ...
                                size(timeSeries, 1), ...
                                size(timeSeries, 2), motifLengths(j), 5,...
                                motifLengths(j), 0);
                        end
                end
                
                % add the motifs to the motif index
                motifIndex{count, 1} = fileList{i};
                motifIndex{count, 2} = motifStart1;
                motifIndex{count, 3} = motifStart2;
                motifIndex{count, 4} = motifLengths(j);
                motifIndex{count, 5} = minDist;
                motifIndex{count, 6} = samplingInterval;
                
                % increment count only if motif data successfully entered
                count = count + 1;
                
            catch ME1
                % if evaluation fails, display file name, motif length, and
                % error message.
                fprintf(['File ' fileList{i} ...
                    '\ngave the following error (with motif length ' ...
                    num2str(motifLengths(j)) '):\n'...
                    ME1.message '\n\n']);
                errorCount = errorCount + 1;
            end
        end
        
    catch ME2
        % if evaluation fails, display file name and error message.
        fprintf(['File ' fileList{i} ...
            '\ngave the following error:\n'...
            ME2.message '\n\n']);
        errorCount = errorCount + 1;
    end
end