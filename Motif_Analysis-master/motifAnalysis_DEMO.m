% This script runs essentially the full set of motif analyses
% starting from projected eigenworm amplitudes to motif discovery and the
% calculation of a motif distance-based feature matrix and affinity
% propagation clustering of mutants.  The analysis and results are
% described in the following paper:
% 
% André E.X. Brown, Eviatar I. Yemini, Laura Grundy, Tadas Jucikas, William
% R. Schafer (2013) A dictionary of behavioral motifs reveals clusters of 
% genes affecting Caenorhabditis elegans locomotion. Proc. Nat. Acad. Sci.
% USA 110:791.  http://dx.doi.org/10.1073/pnas.1211447110
% 
% If you use this script or the associated functions, please cite our
% paper!

% Copyright Medical Research Council 2013
% André Brown, abrown@mrc-lmb.cam.ac.uk, aexbrown@gmail.com
% Released under the terms of the GNU General Public License, see
% readme.txt

 
% To run, add the Motif_Analysis directory to Matlab's path.  Some of the
% associated functions use several mex files to increase speed.  In all 
% cases the c or c++ source is included so if these files don't run on your
% system, you can compile them.  See Matlab's mex files guide for details:
% http://www.mathworks.co.uk/support/tech-notes/1600/1605.html (jan 2013)

% set the root directory where the projected amplitudes are
directory = '/Users/abrown/Andre/wormVideos/Motif_Analysis_SampleData/sample_data/';

% find feature data files and their corresponding worm names.  This
% step relies on a regular expression to parse the file path to get the
% worm name.  This will only work if the projected amps files are arranged
% as directory/mutant name/(optionally other folders)/L, R, or
% unknown (indicates worm side)/projectedAmpsFileName
[fileList, ~] = dirSearch(directory, 'features.mat');

% load the eigenworm projected amplitude data from each feature file,
% interpolate of the NaN values, and save an interpolated time series file
% for each feature file.
minPoints = 13000;
interpNaN(fileList, 'projectedAmpsNoNaN.mat', minPoints, 1, 1)

% re-create fileList looking for interpolated NaN files
[fileList, wormNames] = dirSearch(directory, 'projectedAmpsNoNaN.mat');

% define the parameters for the motif search
numEigWorms = 4;
samplingInterval = 4;
motifLengths = [10 50 150];

% find motif positions for every file in fileList
motifIndex = findMotifsFileList(fileList, numEigWorms, ...
    samplingInterval, motifLengths, 'adaptive', 1);

% keep only the best matching 80% of motifs of each length
[reducedIndex, keptIndicesMatch] = ...
    findBestMatchingMotifs(motifIndex, 'best each length', 0.8);

% using the file paths, motif starts, and sampling interval in the reduced
% index, get the actual full-resolution motifs from the parent files
motifDictionary = index2dict(reducedIndex, numEigWorms);

% get a variety of motifs that represent the full range of motion
% represented in the motifDictionary.
[reducedDictionary, keptIndicesVariable] = ...
    findVariableMotifs(motifDictionary, 'each length', 0.5);

% % plot the motifs in the dictionary (press any key to go to next plot)
% motifPlot(reducedDictionary, numEigWorms, ...
%     max(motifLengths)*samplingInterval, 20, 0, 1)

% use the motif dictionary to create a feature matrix of distances (and
% frequencies if calcFreq == 1)
calcFreq = 0;
featureMatrix = motifFeatureMatrix(fileList, reducedDictionary, ...
    calcFreq, numEigWorms, samplingInterval, 15000, 1);

% sometimes a file doesn't work (e.g. is too short) leading to NaNs in
% feature matrix.  Remove these rows now if present.
nanRows = find(isnan(featureMatrix(:,1)));
featureMatrix(nanRows,:) = [];
wormNames(nanRows) = [];

% normalise the feature matrix (subtract the mean from each column and
% divide by the standard deviation)
meanMatrix = repmat(mean(featureMatrix), size(featureMatrix, 1), 1);
stdMatrix = repmat(std(featureMatrix), size(featureMatrix, 1), 1);
featureMatrixNorm = (featureMatrix - meanMatrix) ./ stdMatrix;

% convert worm names to integers indicating grouped mutant classes
[~, ~, classIntegers] = unique(wormNames);

% use mutual information based criterion to select the least redundant,
% most relevant features (feature indices is a ranking of the features)
numFeatures = size(featureMatrixNorm, 2);
featureIndices = mRMR_feature_select(featureMatrixNorm, classIntegers, ...
    numFeatures);

% if any features are perfectly constant for every sample in a class, the
% covariance matrix will be singular and the Mahalanobis distance can't be
% calculated.  Find any such features and remove them from featureIndices
constantFeatures = findConstantFeatures(featureMatrixNorm, wormNames);

% remove the constant features from featureIndices.  Beware of using
% setdiff for this because it sorts its output!!
featureIndicesNoConst = featureIndices;
for ii = 1:length(constantFeatures)
    indToRemove = find(featureIndicesNoConst == constantFeatures(ii));
    featureIndicesNoConst(indToRemove) = [];
end

% calculate the distance between each of the mutant classes, in this case
% keeping the top 40 features
[distList, uniqueClassNames] = ...
    calcDistList(featureMatrixNorm, wormNames, ...
    featureIndicesNoConst(1:40), 1, 'mahalanobis');

% convert the distance list to a matrix for plotting
N = numel(uniqueClassNames);
distMat = zeros(N, N);
for i = 1:size(distList,1)
    distMat(distList(i, 1), distList(i, 2)) = distList(i, 3);
    distMat(distList(i, 2), distList(i, 1)) = distList(i, 3);
end

% plot the Mahalanobis distance matrix
figure; imagesc(distMat);

% use affinity propagation to cluster mutant classes based on the inverse
% of the distance as the similarity.  The preference for each
% mutant to be an examplar is simply set as the median similarity, as
% suggested by Frey & Dueck.
preference = median(1./distList(:,3));
clusterList = affinityPropagation(distList, uniqueClassNames, preference);

% use bootstrapped feature matrices to determine how often mutant classes
% co-cluster with resampling
numBootstraps = 100;
preferenceFactor = 1;
[clusterFrequency, ~] = ...
    bootstrapAP(featureMatrixNorm, wormNames, ...
    featureIndicesNoConst(1:40), preferenceFactor, numBootstraps, 1, ...
    'mahalanobis');

% keep the top 60% of connections
sortedFreq = sort(cell2mat(clusterFrequency(:, 3)));
connectionThreshold = sortedFreq(round(0.4 * numel(sortedFreq)));

robustConnectionRows = cellfun(@(x) x > connectionThreshold, ...
    clusterFrequency(:, 3));
robustConnections = clusterFrequency(robustConnectionRows, :);

% if there are any mutant classes that don't participate in any robust
% clusters, add them as self-connections with probability 1 so that they
% are still displayed in a newtwork plot by, e.g. Cytoscape

% find the classes that are present in at least one robust cluster
presentClasses = ...
    unique([robustConnections(:, 1); robustConnections(:, 2)]);
% find any missing classes
missingClasses = setdiff(uniqueClassNames, presentClasses);

% make a list of self-connections for any missing classes
selfConnections = cell(numel(missingClasses), 3);
for i = 1:numel(missingClasses)
    selfConnections{i, 1} = missingClasses{i};
    selfConnections{i, 2} = missingClasses{i};
    selfConnections{i, 3} = 0;
end

% add the self connections to the robustConnections network
robustConnections = [robustConnections; selfConnections];

% add the inverse distance information for each robust connection
clusterNetwork = cell(size(robustConnections, 1), 4);

for ii = 1:size(robustConnections, 1)
    % convert the names to indices in unique class names
    classInd1 = find(strcmp(uniqueClassNames, robustConnections{ii, 1})...
        == 1);
    classInd2 = find(strcmp(uniqueClassNames, robustConnections{ii, 2})...
        == 1);
    
    % get the inverse distance from distMat
    invDist = 1/distMat(classInd1, classInd2);
    
    % add to current row of clusterNetwork
    clusterNetwork{ii, 1} = robustConnections{ii, 1};
    clusterNetwork{ii, 2} = robustConnections{ii, 2};
    clusterNetwork{ii, 3} = robustConnections{ii, 3};
    clusterNetwork{ii, 4} = invDist;
end