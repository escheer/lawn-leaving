function constantFeatures = findConstantFeatures(featureMatrix, classNames)

% FINDCONSTANTFEATURES Identify any features that are constant across all
% samples in a class.
% 
% Input
%   featureMatrix  - a feature matrix with samples along rows and
%                    features down columns.
%   classNames     - a list of class names for each row of featureMatrix
% 
% Output
% constantFeatures - a vector containing the feature (column) indices in
%                    feature matrix for which the feature values are
%                    constant for every row in at least one class.
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
constantFeatures = [];

% loop through classes
for i = 1:numel(uniqueClassNames);
    % get the data for this class name
    name = uniqueClassNames{i};
    rows = strcmp(classNames, name);
    data = featureMatrix(rows, :);
    
    % find any features that are totally constant in data and add them to
    % constantFeatures.  Constant features should be rare, so don't worry 
    % about cost of growing constantFeatures in the loop.
    constantFeatures = [constantFeatures find(all( diff(data) == 0 ))];
end

% remove repetitions
constantFeatures = unique(constantFeatures);