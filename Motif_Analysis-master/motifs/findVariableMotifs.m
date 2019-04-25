function [reducedDictionary keptIndices] = ...
    findVariableMotifs(motifDictionary, selectionCriterion, fractionToKeep)

% FINDVARIABLEMOTIFS Given a motif index, find the motifs with the
% highest variability.  Only the top fractionToKeep are output in the
% reduced index.
%
% Input
%   motifDictionary    - A cell array containing the parent file names, the
%                        actual motifs, and their start locations at the
%                        original time series resolution (i.e. no longer
%                        down sampled)
%   selectionCriterion - string specifying how indices are searched for.
%                        Can be 'overall' or 'each length' to require that
%                        selected motifs come from each possible length.
%   fractionToKeep     - a number from 0 to 1 specifying what fraction of
%                        the input motif list should be kept as best
%                        matches
%
% Output
%   reducedDictionary  - reduced motif dictionary that includes only the
%                        most variable motifs
%   keptIndices        - the row indices in the original motif dictionary
%                        that were kept
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


% initialise a vector for the mean standard deviation of the motifs
meanSTD = NaN(size(motifDictionary, 1), 1);

% loop through all motifs in the dictionary and find their standard
% deviation
for ii = 1:size(motifDictionary, 1)
    % get the mean standard deviation of the two motifs
    motif1 = motifDictionary{ii, 3};
    motif2 = motifDictionary{ii, 5};
    
    meanSTD(ii) = mean( [std(motif1, 0, 2);  std(motif2, 0, 2)] );
end


% use the appropriate selection case
switch selectionCriterion
    
    % simply take the best matching cases in the motif index
    case 'overall'
        
        % sort based on distance between motifs
        [~, sortInds] = sort(meanSTD);
        
        % sample uniformly across the sorted indices so that the output
        % has a variety of motifs with different amounts of motion
        stepSize = round(length(sortInds) * fractionToKeep);
        variableInds = sortInds(1:stepSize:end);
        
        % check if there are any rows to keep
        if numRowsToKeep >= 1
            %get the reduced index and reduced file list
            keptIndices = variableInds;
            reducedDictionary = motifDictionary(variableInds, :);
        else
            keptIndices = [];
            reducedDictionary = cell(0,5);
        end
        
        % take the best matches from each length judged separately
    case 'each length'
        
        % get the motif lengths
        motifLengths = NaN(size(motifDictionary, 1), 1);
        for i = 1:size(motifDictionary, 1)
            motifLengths(i) = size(motifDictionary{i, 3}, 2);
        end
        
        % get the number of motif lengths
        [uniqueLengths, ~] = unique(motifLengths);
        
        % initialise output. Rather complicated expression for reduced
        % index length is to account for rounding number of rows to include
        % for each length, which is potentially different from
        % round(size(motifIndex, 1) * fractionToKeep)
        reducedDictionary = ...
            cell(round( size(motifDictionary, 1)/length(uniqueLengths) * ...
            fractionToKeep ) * length(uniqueLengths), ...
            size(motifDictionary, 2));
        keptIndices = ...
            NaN(round( size(motifDictionary, 1)/length(uniqueLengths) * ...
            fractionToKeep ) * length(uniqueLengths), 1);
        
        % process each motif length separately
        for i = 1:length(uniqueLengths)
            % get indices of motifs with the current length
            motifInds = find(motifLengths == uniqueLengths(i));
            
            % sort based on distance between motifs
            [~, sortInds] = sort(meanSTD(motifInds));
            
            % sample uniformly across the sorted indices so that the output
            % has a variety of motifs with different amounts of motion
            stepSize = round(1/fractionToKeep);
            variableInds = motifInds(sortInds(1:stepSize:end));
            
            % get the reduced index and file list for this length
            reducedDictionary(1 + (i - 1)*length(variableInds) : ...
                i*length(variableInds), :) = motifDictionary(variableInds, :);
            keptIndices(1 + (i - 1)*length(variableInds) : ...
                i*length(variableInds), :) = variableInds;
        end
        
end