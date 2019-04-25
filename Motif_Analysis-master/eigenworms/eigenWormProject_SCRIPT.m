% load angle arrays, flip for dorsoventral consistency, project onto
% eigenworms, and interpolate over NaN values, dropping the NaN values at
% the beginning and end to avoid unrealistic extrapolations

% load the eigenworms
load('/Users/abrown/Andre/wormVideos/results-12-05-10/eigenwormsN2.mat')

% what's the minimum number of non-NaN frames that are allowed?
minPoints = 10000;

% set the root directory
directory = ['/Users/abrown/Andre/wormVideos/results-12-05-10/'...
    'wild-isolates/'];

% get the list of files
[fileList, wormNames] = dirSearch(directory, 'angleArray.mat');

% loop through files
for ii = 1:numel(fileList)
    disp(ii/numel(fileList))
    
    % load the current file
    angles = cell2mat(struct2cell( load(fileList{ii}) ));
    
    % if the worm is in the "left" folder, invert all the angles.
    if ~isempty(strfind(fileList{ii}, '/L/'))
        angles = angles * -1;
    end
    
    % get the projected amplitudes
    numEigWorms = 20;
    timeSeries = eigenWormProject(eigenWorms, angles, numEigWorms)';
    
    % remove leading and trailing NaN values
    % is the first point NaN?
    if isnan(timeSeries(1, 1))
        % get the end of the starting NaN segment
        nanEnd = find(~isnan(timeSeries(1, :)), 1, 'first') - 1;
        
        % drop these values
        timeSeries(:, 1:nanEnd) = [];
    end
    
    % is the last point NaN?
    if isnan(timeSeries(1, end))
        % get the start of the final NaN segment
        nanStart = find(~isnan(timeSeries(1, :)), 1, 'last') + 1;
        
        % drop these values
        timeSeries(:, nanStart:end) = [];
    end
    
    % check that sufficient frames are not NaN. Only check first
    % dimension since if any value is NaN in a frame all dimensions should
    % be.
    if sum(~isnan(timeSeries(1, :))) >= minPoints
        
        % initialise timeSeriesNoNaN
        projectedAmpsNoNaN = timeSeries;
        
        % interpolate over NaN values
        for j = 1:size(timeSeries, 1)
            pAmp = timeSeries(j, :);
            pAmp(isnan(pAmp)) = interp1(find(~isnan(pAmp)),...
                pAmp(~isnan(pAmp)), find(isnan(pAmp)),'linear', 'extrap');
            projectedAmpsNoNaN(j, :) = pAmp;
        end
        
        % save interpolated data
        save(strrep(fileList{ii}, 'angleArray', 'projectedAmpsNoNaN'), ...
            'projectedAmpsNoNaN')
    end
end