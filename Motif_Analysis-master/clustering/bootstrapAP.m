function [clusterProbabilities clusterSizes] = ...
    bootstrapAP(featureMatrix, classNames, featureInds, ...
    preferenceFactor, numberOfBootstraps, display, distType)

% BOOTSTRAPAP Use resampling with replacement on featureMatrix and
% repeat affinity propagation clustering for each resampled data set.  
% Record fraction of times that each mutant class clusters with each other 
% class.
%
% Input
%   featureMatrix        - a feature matrix with samples along rows and
%                          features down columns.
%   classNames           - a list of class names for each row of
%                          featureMatrix
%   featureInds          - a vector of indices indicating which features
%                          (columns) of featureMatrix should be used in
%                          calculating the distance
%   preferenceFactor     - Preference is set to preferenceFactor/
%                          median(similarity).  From the documentation for 
%                          apcluster:
%                          "a real-valued N-vector. p(i) indicates the
%                          preference that data point i be chosen as an
%                          exemplar. Often a good choice is to set all
%                          preferences to median(s); the number of clusters
%                          identified can be adjusted by changing this
%                          value accordingly. If 'p' is a scalar, APCLUSTER
%                          assumes all preferences are that shared value."
%   numberOfBootstraps   - how many times should featureMatrix be resampled
%                          and clustered?
%   display              - logical.  Should the progress of the calculation
%                          be displayed in the command window?
%   distType             - string indicating whether distance should be
%                          'euclidean' or 'mahalanobis'
%
% Output
%   clusterProbabilities - for each class name combination, what fraction
%                          of the bootstrap runs resulted in co-clustering?
%                          (e.g. 0 indicates this combination never
%                          clustered together, 0.5 indicates that this pair
%                          clustered together in half of the runs)
%   clusterSizes         - a cell array of vectors, each of which contains
%                          the sizes of the clusters created during each
%                          bootstrap round.  numel(clusterSizes{i}) gives
%                          the cluster number for the ith round.
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


% check arguments
if preferenceFactor <= 0
    error('preferenceFactor must be greater than 0.')
end

% initialise clusterNumber
clusterSizes = cell(numberOfBootstraps, 1);

% get the unique class names
uniqueClassNames = unique(classNames);
N = numel(uniqueClassNames);

% initialise clusterProbabilities
clusterProbabilities = cell((N^2 - N)/2, 3);

% fill in class names in clusterProbabilities
count = 1;
for i = 1:N
    for j = i + 1:N
        clusterProbabilities(count, :) = ...
            {uniqueClassNames{i}, uniqueClassNames{j}, 0};
        count = count + 1;
    end
end

% if number of bootstraps is 0, do not resample
if numberOfBootstraps == 0
    resampleFeatMat = 0;
    % add 1 to number of bootstraps so that one iteration of clustering is
    % done below
    numberOfBootstraps = 1;
end

% resample and re-cluster featureMatrix specified number of times
for i = 1:numberOfBootstraps
    % display progress in command window if requested
    if display
        disp(i/numberOfBootstraps)
    end
    
    % resample from the feature matrix
    featMatResamp = resampleFeatureMatrix(featureMatrix, classNames);
    
    % calculate Mahalanobis distance using resampled data
    [distList, ~] = ...
        calcDistList(featMatResamp, classNames, featureInds, 0, distType);
    
    % distList should normally not contain any NaN values, but if the
    % pooled covariance matrix for a given comparison is badly conditioned
    % even after removing repeating features, there could be some NaNs.
    % This could happen if you are comparing two mutant classes each with
    % two or fewer samples. In case this happens, remove rows with NaNs.
    distList(isnan(distList(:,3)), :) = [];
    
    % use affinity propagation to cluster mutant classes based on the
    % inverse of the Mahalanobis distance as the similarity.  The
    % preference for each mutant to be an examplar is simply set as the
    % median similarity by default, as suggested by Frey & Dueck, but this
    % can be modified by preferenceFactor if more or fewer clusters are
    % desired.
    preference = median(1./distList(:,3)) * preferenceFactor;
    clusterList = ...
        affinityPropagation(distList, uniqueClassNames, preference);
    
    % get the names of the exemplars determined by the clustering algorithm
    exemplars = unique(clusterList(:, 2));
    
    % initialise a vector to hold the cluster sizes for the current
    % bootstrap run
    currentClusterSizes = [];
    
    % determine connections in cluster list and add them to
    % clusterProbabilities
    for j = 1:numel(exemplars)
        % find the current exmplar and the corresponding rows in cluster
        % list
        currentExemplar = exemplars{j};
        exemplarInds = strcmp(clusterList(:,2), currentExemplar);
        
        % get the names of each class in the current cluster
        clusterMembers = clusterList(exemplarInds, 1);
        
        % get the distribution of cluster sizes
        currentClusterSizes = [currentClusterSizes numel(clusterMembers)];
        
        %         % remove the exemplar from the list of cluster members (we don't
        %         % count a class clustering with itself)
        %         selfMatchInd = strcmp(clusterMembers, currentExemplar);
        %         clusterMembers(selfMatchInd) = [];
        
        % add each of the associations represented by the current cluster
        % to clusterProbabilities.  If the clusterMembers consisted only of
        % the exmemplar itself, clusterMembers should now be empty and this
        % step will simply be skipped.
        for k = 1:numel(clusterMembers)
            
            % find the current cluster member in clusterProbabilities
            currentMemberRows = ...
                find(strcmp(clusterProbabilities(:,1), clusterMembers{k})...
                == 1);
            
            % add all the co-clustering classes.  Go from k + 1 to end to
            % avoid counting self associations or over counting proper
            % associations
            for m = k + 1:numel(clusterMembers)
                % get the mth co-clustering class
                coClusterer = clusterMembers{m};
                
                % find the row corresponding to the current association
                coClustererRow = ...
                    strcmp(clusterProbabilities(currentMemberRows, 2),...
                    coClusterer);
                
                % make sure that there is no more than one matching row
                if sum(coClustererRow) > 1
                    error(['There is more than one matching entry in '...
                        'clusterProbabilities. Is it poorly initialised?'])
                end
                
                % find the index in clusterProbabilities of the matching
                % row
                matchingRow = currentMemberRows(coClustererRow);
                
                % increment the number of matches
                clusterProbabilities{matchingRow, 3} = ...
                    clusterProbabilities{matchingRow, 3} + 1;
            end
        end
    end
    % add the cluster size vector to clusterSizes
    clusterSizes{i} = currentClusterSizes;
end


% divide number of co-clusterings by number of bootstraps.
for i = 1:size(clusterProbabilities, 1)
    clusterProbabilities{i, 3} = ...
        clusterProbabilities{i, 3} / numberOfBootstraps;
end
