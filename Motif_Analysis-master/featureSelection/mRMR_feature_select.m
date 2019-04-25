function featureIndices = ...
    mRMR_feature_select(featureMatrix, classIntegers, numFeatures)

% MRMR_FEATURE_SELECT Discretise feature matrix using SAX representation of
% Lin and Keogh then uses Peng's minimum redundancy, maximum relevance to
% select features from the feature matrix.
%
% References:
%
%   Lin, J., Keogh, E., Lonardi, S. & Chiu, B.
%   "A Symbolic Representation of Time Series, with Implications for
%   Streaming Algorithms."
%   In proceedings of the 8th ACM SIGMOD Workshop on Research Issues in
%   Data Mining and Knowledge Discovery. San Diego, CA. June 13, 2003.
%
%   Hanchuan Peng, Fuhui Long, and Chris Ding
%   "Feature selection based on mutual information: criteria of
%   max-dependency, max-relevance, and min-redundancy,"
%   IEEE Transactions on Pattern Analysis and Machine Intelligence,
%   Vol. 27, No. 8, pp.1226-1238, 2005.
%
%
% Input
%   featureMatrix  - a feature matrix with samples along rows and features
%                    down columns.
%   classIntegers  - a series of integers to indicate which class the rows
%                    in featureMatrix belong to. 1 for first class, 2 for
%                    second class, etc.
%   numFeatures    - The number of features that should be included in the
%                    ranked list 
%
% Output
%   featureIndices - feature indices ranked by mRMR.  If using e.g. 5
%                    features for classification, use
%                    featureMatrix(:, featureIndices(1:5))
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


% get size of feature matrix
[nRows nCols] = size(featureMatrix);

% set number of levels to use in discretisation
nLevels = 10;

discreteFeatureMatrix = NaN(nRows, nCols);

for i = 1:nCols
    % calculate symbolic approximation of data
    [symbolic_data, ~] =  ...
        timeseries2symbol(featureMatrix(:, i),nRows, nRows, nLevels, 1);
    discreteFeatureMatrix(:,i) = symbolic_data - ceil(nLevels/2);
end

% use mRMR to rank features
featureIndices = ...
    mrmr_miq_d(discreteFeatureMatrix, classIntegers, numFeatures);