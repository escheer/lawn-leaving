function [reducedIndex keptIndices] = ...
    findBestMatchingMotifs(motifIndex, selectionCriterion, fractionToKeep)

% FINDBESTMATCHINGMOTIFS Given a motif index, find the best matching motifs
% according to the specified selection criterion (best overall, best each
% length, best each class).  Only the top fractionToKeep are output in the
% reduced index.
% 
% Input
%   motifIndex         - an index of motifs (includes starts of motifs, 
%                        their distance from each other, and their length)
%   selectionCriterion - string specifying how best matches are defined.
%                        Can be 'best overall' or 'best each length'
%   fractionToKeep     - a number from 0 to 1 specifying what fraction of
%                        the input motif list should be kept as best
%                        matches
% 
% Output
%   reducedIndex       - the reduced motif index that includes only the
%                        best matches
%   keptIndices        - the row indices in the original motif index that
%                        were kept (potentially useful for reducing other
%                        things like a class name list)
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


% use the appropriate selection case
switch selectionCriterion
    
    
    % simply take the best matching cases in the motif index
    case 'best overall'
        
        % sort based on distance between motifs
        [~, sortInds] = sort(cell2mat(motifIndex(:, 5)));
        
        % get number of rows to keep
        numRowsToKeep = round(size(motifIndex, 1) * fractionToKeep);
        
        % check if there are any rows to keep
        if numRowsToKeep >= 1
            %get the reduced index and reduced file list
            keptIndices = sortInds(1:numRowsToKeep);
            reducedIndex = motifIndex(keptIndices, :);
        else
            reducedIndex = cell(0,6);
            keptIndices = [];
        end
        
    % take the best matches from each length judged separately
    case 'best each length'
        
        % get the number of motif lengths
        [motifLengths, ~] = unique(cell2mat(motifIndex(:, 4)));
        
        % initialise output. Rather complicated expression for reduced
        % index length is to account for rounding number of rows to include
        % for each length, which is potentially different from
        % round(size(motifIndex, 1) * fractionToKeep)
        reducedIndex = ...
            cell(round( size(motifIndex, 1)/length(motifLengths) * ...
            fractionToKeep ) * length(motifLengths), ...
            size(motifIndex, 2));
        keptIndices = ...
            NaN(round( size(motifIndex, 1)/length(motifLengths) * ...
            fractionToKeep ) * length(motifLengths), 1);
        
        % process each motif length separately
        for i = 1:length(motifLengths)
            % get indices of motifs with the current length
            motifInds = find(cell2mat(motifIndex(:, 4)) == ...
                motifLengths(i));
            
            % sort based on distance between motifs
            [~, sortInds] = sort(cell2mat(motifIndex(motifInds, 5)));
            
            % take rankedInds and find corresponding indices in original
            % motifLibrary (which includes motifs of all lengths, not just
            % the currently searched one)
            bestMatchInds = ...
                motifInds(sortInds(1:round(length(sortInds) * ...
                fractionToKeep)));
            
            % get the reduced index and file list for this length
            reducedIndex(1 + (i - 1)*length(bestMatchInds) : ...
                i*length(bestMatchInds), :) = motifIndex(bestMatchInds, :);
            keptIndices(1 + (i - 1)*length(bestMatchInds) : ...
                i*length(bestMatchInds), :) = bestMatchInds;
        end
end