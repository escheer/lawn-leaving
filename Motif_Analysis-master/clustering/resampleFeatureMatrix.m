function featMatResamp = resampleFeatureMatrix(featureMatrix, classNames)

% RESAMPLEFEATUREMATRIX Take a feature matrix and resample with replacement
% from each of the classes to make a new bootstrapped feature matrix.
% 
% Input
%   featureMatrix  - a feature matrix with samples along rows and features
%                    down columns.
%   classNames     - a list of class names for each row of featureMatrix
% 
% Output
%   featMatResamp  - the resampled feature matrix
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


% get each class integer (this approach allows for things like 
% non-consecutive class integers in case that arises)
uniqueClassNames = unique(classNames);

% generate list of row indices by resampling separately for each class and
% combining the results
resampledRows = NaN(size(featureMatrix, 1), 1);
for i = 1:length(uniqueClassNames)
    % find the original rows in featureMatrix of the current class
    currentRows = find(strcmp(classNames, uniqueClassNames{i}) == 1);
    
    % resample with replacement from currentRows and add the resampled data
    % to resampledRowInds
    resampledRows(currentRows) = ...
        randi([currentRows(1), currentRows(end)], 1, length(currentRows));
end

% use the resampled rows to make a resampled version of featureMatrix
featMatResamp = featureMatrix(resampledRows, :);    