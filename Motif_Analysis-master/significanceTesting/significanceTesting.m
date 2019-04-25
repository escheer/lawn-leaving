function [t2Stats pVals] = ...
    significanceTesting(featureMatrix, classNames, featureInds, refName,...
    equalVariance, permutations, display)

% SIGNIFICANCETESTING significance testing on a feature matrix using the
% Hotelling T^2 statistic (calculated with a shrinkage estimate of the
% covariance matrix) and a permutation test for determining the p-value for
% the comparison of data with refName against all other data classes in the
% feature matrix.
%
% Input
%   featureMatrix - a feature matrix with samples along rows and features 
%                   down columns.
%   classNames    - a list of class names for each row of featureMatrix
%   featureInds   - a vector of indices indicating which features (columns)
%                   of featureMatrix should be used in calculating the 
%                   Mahalanobis distance
%   refName       - the reference name that all other classes should be 
%                   compared against
%   equalVariance - Logical.  Are the variances of the two data sets equal?
%   permutations  - The number of permutations to perform in the
%                   permutation test.  Set to 0 to skip permutation test.
%                   In this case pVal will be returned as NaN.
%   display       - Logical.  Should the progress be displayed in the
%                   command window?
%
% Output
%   t2Stats       - A vector of Hotelling T^2 statistics for each 
%                   comparison
%   pVals         - The p-values for each comparison derived from the
%                   permutation test of the T^2 statistics
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

% get the unique class names
uniqueClassNames = unique(classNames);

% initialisation
t2Stats = NaN(numel(uniqueClassNames), 1);
pVals = NaN(numel(uniqueClassNames), 1);

% get the rows in featureMatrix corresponding to refName
refRows =  strcmp(classNames, refName);
refData = featureMatrix(refRows, featureInds);

% loop through uniqueClassNames and compare them to the reference data
for i = 1:numel(uniqueClassNames)
    % if desired, display the progress in the command window
    if display
        disp(i/numel(uniqueClassNames))
    end
    
    testName = uniqueClassNames{i};
    
    % don't bother testing the reference name against itself
    if ~strcmp(testName, refName)
        % get the test data set
        testRows = strcmp(classNames, testName);
        testData = featureMatrix(testRows, featureInds);
        
        [t2Stats(i) pVals(i)] = ...
            hotellingT2(refData, testData, equalVariance, permutations);
    end
end