function interpNaN(fileList, outputName, minPoints, trim, display)

% INTERPNAN Load each time series file from fileList and interpolate over
% the NaN values using linear interpolation and extrapolation.
% Extrapolation is needed in case the first or last values are NaN. If
% every point in a file is NaN, then interpNaN skips the file. The
% interpolated files are saved to the same location as the un-interpolated
% feature file and called outputFilename
%
% Input
%   fileList   - List of the full paths of the feature files to be
%                processed.  Note that the projected amplitudes data is
%                only part of this file.
%   outputName - The file name to be used for the interpolated data.  Will
%                be placed in the same directory as the file in fileList.
%   minPoints  - Files with fewer than minPoints entries that are non NaN
%                will be skipped.
%   trim       - Logical.  Should leading and trailing NaN values be
%                dropped before interpolating?
%   display    - Logical.  Should the progress be displayed in the
%                command window?
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

% loop through fileList
for i = 1:size(fileList, 1)
    % display progress
    if display == 1
        disp( ['interpNaN progress: ' num2str(i/size(fileList, 1))] )
    end
    
    % get the current fileName path
    filePath = fileList{i};
    
    % load the time series to examine
    worm = cell2mat(struct2cell(load( filePath, 'worm' )));
    timeSeries = worm.posture.eigenProjection;
    
    % optionally, remove leading and trailing NaN values (these are
    % "unanchored" and if they are long can lead to extreme outlying
    % values)
    if trim == 1
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
    end
    
    % check that at least two values are not NaN. Only check first
    % dimension since if any value is NaN in a frame all dimensions should
    % be.
    if sum(~isnan(timeSeries(1, :))) >= minPoints
        
        % initialise timeSeriesNoNaN
        timeSeriesNoNaN = timeSeries;
        
        % interpolate over NaN values
        for j = 1:size(timeSeries, 1)
            pAmp = timeSeries(j, :);
            pAmp(isnan(pAmp)) = interp1(find(~isnan(pAmp)),...
                pAmp(~isnan(pAmp)), find(isnan(pAmp)),'linear', 'extrap');
            timeSeriesNoNaN(j, :) = pAmp;
        end
        
        % save interpolated data
        
        % get the directory from fileList and save to output file. Slashes
        % must be different for PC vs. mac OS, so check if you're on a pc.
        if ispc
            directory = char(regexpi(filePath, '.+(?=\\)', 'match'));
            save([directory '\' outputName], 'timeSeriesNoNaN')
        else
            directory = char(regexpi(filePath, '.+(?=/)', 'match'));
            save([directory '/' outputName], 'timeSeriesNoNaN')
        end
    end
end