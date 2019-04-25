function motifDictionary = index2dict(motifIndex, numDimensions)

% INDEX2DICT  Takes a motifIndex cell array and returns the corresponding
% motif dictionary.  It simply goes into the files listed in each row of
% the motifIndex and uses the start locations, sampling interval, and motif
% length to extract the actual motifs.
%
% Input
%   motifIndex      - A motif index cell array with parent file names,
%                     motif start indices, lengths, and sampling intervals.
%   numDimensions   - The number of eigen worms to consider in motif
%                     finding
%   motifDictionary - A cell array containing the parent file names, the
%                     actual motifs, and their start locations at the
%                     original time series resolution (i.e. no longer down
%                     sampled)
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


% initialise motifDictionary
motifDictionary = cell(size(motifIndex, 1), 5);

% loop through rows of motif index
for i = 1:size(motifIndex, 1)
    
    % check for empty entries in motifIndex
    if ~isempty(motifIndex{i, 1})
        
        % load the motif parent time series
        timeSeries = cell2mat(struct2cell(load( motifIndex{i, 1} )));
        
        % get info from motifIndex
        start1 = motifIndex{i, 2};
        start2 = motifIndex{i, 3};
        motifLength = motifIndex{i, 4};
        samplingInterval = motifIndex{i, 6};
        
        % extract the two motifs from the parent
        motif1 = timeSeries( 1:numDimensions,...
            (start1 - 1) * samplingInterval + 1:...
            (start1 + motifLength - 2) * samplingInterval + 1);
        motif2 = timeSeries( 1:numDimensions,...
            (start2 - 1) * samplingInterval + 1:...
            (start2 + motifLength - 2) * samplingInterval + 1);
        
        % add the motifs to the motif dictionary
        motifDictionary{i, 1} = motifIndex{i, 1};
        motifDictionary{i, 2} = start1*samplingInterval;
        motifDictionary{i, 3} = motif1;
        motifDictionary{i, 4} = start2*samplingInterval;
        motifDictionary{i, 5} = motif2;
    end
end