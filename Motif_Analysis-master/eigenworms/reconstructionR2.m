function rSquared = reconstructionR2(fileList, eigenWorms, ...
    numDimensions, display)

% RECONSTRUCTIONR2  Use eigenworms to reconstruct angle data and skeleton
% data to compare with the original.  Calculate correlation between
% original and reconstructed data as a quality check.
%
% Input
%   fileList      - The full path of projected amplitude files to analyse.
%                   The projected amplitude files must be called
%                   projectedAmpsNoNaN.mat.  The corresponding skeleton 
%                   angles files must be in the same directory with 
%                   filename eigenAngles.mat
%   eigenWorms    - the basis eigenWorms that were used in the
%                   decomposition
%   numDimensions - The number of eigen worms to consider in motif
%                   finding
%   display       - Logical.  Should the progress be displayed in the
%                   command window?
%
% Output
%   rSquared      - A vector of R^2 values for the correlation between the
%                   raw skeleton angles and the reconstructed angles from
%                   the eigenworm decomposition.  Values are in the same
%                   order as specified by fileList.
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


% initialisation
rSquared = NaN(size(fileList, 1), 1);

% loop through all files in the file list
for i = 1:size(fileList, 1)
    
    % display the progress, if requested
    if display == 1
        disp(i/size(fileList, 1))
    end
    
    % load the time series
    timeSeries = cell2mat(struct2cell(load( fileList{i} )))';
    
    % load the angle array
    angleFileName = strrep(fileList{i}, 'projectedAmpsNoNaN.mat', ...
        'angleArray.mat');
    angleArray = cell2mat(struct2cell(load( angleFileName )));
    
    
    reconstructedAngle = zeros(size(angleArray));
    %reconstruct worm skeleton angles (add up eigenworms times amplitudes)
    for j = 1:size(angleArray, 1)
        amplitudes = timeSeries(j, :);
        for k = 1:numDimensions
            a = amplitudes(k) .* eigenWorms(k, :); %amplitude times corresponding eigenWorm
            reconstructedAngle(j, :) = reconstructedAngle(j, :) + a;
        end
    end
    
    %calculate correlation between angleArray and
    %reconstructedAngleArray to check for degree of agreement
    anglesNoNan = angleArray(~isnan(angleArray));
    reconAngles = reconstructedAngle(~isnan(angleArray));
    R = corrcoef(anglesNoNan,reconAngles);
    rSquared(i) = R(1,2)^2;
end

%make plot of rSquared values to check for possible errors
plot(rSquared)
