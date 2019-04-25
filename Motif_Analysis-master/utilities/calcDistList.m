function [distList uniqueClassNames] = ...
    calcDistList(featureMatrix, classNames, featureInds, display, distType)

% CALCDISTLIST Converts a data matrix with samples = rows and features =
% columns to a distance list.  Each row of the distance list gives the
% "source" and "target" indices and their distance.
%
% Input
%   featureMatrix    - a feature matrix with samples along rows and
%                      features down columns.
%   classNames       - a list of class names for each row of featureMatrix
%   featureInds      - a vector of indices indicating which features
%                      (columns) of featureMatrix should be used in
%                      calculating the distance
%   display          - logical.  Should the progress of the calculation be
%                      displayed in the command window?
%   distType         - string indicating whether distance should be
%                      'euclidean' or 'mahalanobis'
%
% Output
%   distList         - a matrix in which row gives the "source" and
%                      "target" indices and their distance.
%                      Indices rather than classNames are used because this
%                      is the format needed for apcluster.
%   uniqueClassNames - a list of class names that correspond to the indices
%                      in distList, so that if the source index is i in
%                      distList then the corresponding class name is
%                      uniqueClassNames{i}
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


% check type argument
if ~strcmp(distType, 'euclidean') && ~strcmp(distType, 'mahalanobis')
    error('The distType must be either euclidean or mahalanobis.')
end

% get the unique class names
uniqueClassNames = unique(classNames);
N = numel(uniqueClassNames);

% initialise distance list
distList = NaN((N^2 - N)/2, 3);

% pre-compute the covariance matrix and mean of each data set to speed up
% calculation of Mahalanobis distance
if strcmp(distType, 'mahalanobis')
    % initialise
    covFullSet = NaN(N, length(featureInds), length(featureInds));
    meanFullSet = NaN(N, length(featureInds));
    
    for ii = 1:N
        % get the current data set
        name = uniqueClassNames{ii};
        rows = strcmp(classNames, name);
        data = featureMatrix(rows, featureInds);
        
        % calculate the mean and covariance matrix
        meanFullSet(ii, :) = mean(data);
        
        [S, ~, ~] = covshrinkKPM(data, 0);
        covFullSet(ii, :, :) = S;
    end
end

% index for counting through distList
count = 1;
for i = 1:N
    % display progress, if requested
    if display
        disp(i/N)
    end
    
    % get the test data set
    name1 = uniqueClassNames{i};
    rows1 = strcmp(classNames, name1);
    data1 = featureMatrix(rows1, featureInds);
    
    % check that there are at least two rows in data1 if the distType is
    % Mahalanobis.  (Even one row is fine for Euclidean distance)
    if (strcmp(distType, 'euclidean') && size(data1, 1) > 0) || ...
            (strcmp(distType, 'mahalanobis') && size(data1, 1) > 1)
        
        % compare the current test data set against all others.  Because of
        % symmetry of distance metric, only compute for i + 1 to N.
        for j = i + 1:N
            name2 = uniqueClassNames{j};
            rows2 = strcmp(classNames, name2);
            data2 = featureMatrix(rows2, featureInds);
            
            % check that there are at least two rows in data1 if the 
            % distType is Mahalanobis.  (Even one row is fine for Euclidean
            % distance)
            if (strcmp(distType, 'euclidean') && size(data1, 1) > 0) || ...
                    (strcmp(distType, 'mahalanobis') && size(data1, 1) > 1)
                if strcmp(distType, 'euclidean')
                    % calculate the Euclidean distance between data1 and
                    % data2
                    distance = groupEuclideanDistance(data1, data2);
                elseif strcmp(distType, 'mahalanobis')
                    % calculate the Mahalanobis distance between data1 and 
                    % data2.
                    distance = ...
                        mahaDistCOV(squeeze(covFullSet(i, :, :)), ...
                        squeeze(covFullSet(j, :, :)), sum(rows1), ...
                        sum(rows2), meanFullSet(i, :), ...
                        meanFullSet(j, :), 1);
                end
                % add indices and distance to distList and increment
                % count.  i and j are added to distList instead of
                % name1 and name2 because this is the format needed for
                % apcluster.
                distList(count, :) = [i, j, distance];
                count = count + 1;
            else
                % increment count to keep alignment between names and
                % distances
                count = count + 1;
            end
        end
    else
        % increment count to keep name alignment with distances
        count = count + length(i + 1:N);
    end
end
