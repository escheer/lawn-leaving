function clusterList = ...
    affinityPropagation(distList, uniqueClassNames, preference)

% AFFINITYPROPAGATION Use affinity propagation to cluster worms using 
% the input distance metric to construct a similarity measure.  Uses the
% apcluster function from Frey and Dueck.
% 
% Input
%   distList         - a matrix in which row gives the "source" and  
%                      "target" indices and their distance.
%   uniqueClassNames - a list of class names for each row of featureMatrix
%   preference       - from the documentation for apcluster:
%                      "a real-valued N-vector. p(i) indicates the 
%                      preference that data point i be chosen as an  
%                      exemplar. Often a good choice is to set all 
%                      preferences to median(s); the number of clusters
%                      identified can be adjusted by changing this value
%                      accordingly. If 'p' is a scalar, APCLUSTER assumes 
%                      all preferences are that shared value."
%
% Output
%   clusterList      - a cell array with two columns.  The first column is
%                      a list of all the unique class names and the second
%                      column is that class's exemplar from the clustering.
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


% convert distance to similarity using inverse
simList = distList;
simList(:,3) = 1./simList(:,3);

N = numel(uniqueClassNames);
% convert the similarity list to a matrix (apcluster accepts similarity
% lists directly, but just converts them to a matrix anyway).  Diagonal
% elements are set to similarity of -realmax (i.e. 1/distance, with
% self-distance being zero)
simMat = -realmax * ones(N, N);
for i = 1:size(simList,1)
    simMat(simList(i, 1), simList(i, 2)) = simList(i, 3);
end

% run affinity propagation function
[exemplarIndices, ~, ~, ~] = apcluster(simMat, preference, 'nonoise', 1);

% sort exemplarIndices for easy viewing
[sortedExemplars sortingOrder] = sort(exemplarIndices);

% make list of worm names (left column) with their corresponding exemplars
% (right column)
clusterList = [uniqueClassNames(sortingOrder) ...
    uniqueClassNames(sortedExemplars)];