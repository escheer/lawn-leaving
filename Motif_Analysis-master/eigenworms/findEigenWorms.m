function [eigenWorms, eigenVals] = findEigenWorms(angleArray, numEigWorms, verbose)

%FINDEIGENWORMS Find the eigenvectors of the angleArray covariance matrix
%
%   [EIGENWORMS, EIGENVALS] = FINDEIGENWORMS(ANGLEARRAY, NUMEIGWORMS, VERBOSE)
%
%   Input:
%       angleArray  - a numFrames by numSkelPoints - 1 array of skeleton
%                     tangent angles that have been rotated to have a mean
%                     angle of zero.  For best result, this angle array
%                     should probably be taken from a couple of hours of
%                     tracking data.  Several 15 minute angleArrays can be
%                     concatenated.
%       numEigWorms - the number of eigenworms to use in the projection.
%                     We simply take the first numEigWorms dimensions from
%                     eigenWorms and project onto those.
%       verbose     - flag indicating whether plots should be generated.
%
%   Output:
%       eigenWorms  - a matrix with the eigenvectors of the covariance
%                     matrix of angleArray ranked according to the amount
%                     variance they account for (most variance is first)
% 
% 
% Copyright Medical Research Council 2013
% Andr� Brown, andre.brown@csc.mrc.ac.uk, aexbrown@gmail.com
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


% check for all-NaN input data
if isempty(~isnan(angleArray))
    error('angleArray contains only NaN values.')
end

% only use non-NaN frames for calculating the covariance matrix
covarianceMat = cov(angleArray,0,'omitrows');

% get the eigenvectors and eigenvalues of the covariance matrix
[M, eVals] = eig(covarianceMat);

% sort the eigenvalues
eVals = sort(diag(eVals),'descend');

% keep the numEigWorms dimensions that capture most of the variance
eigenWorms = M(:, end:-1:end - numEigWorms + 1)';
eigenVals = eVals(1:numEigWorms);

if verbose
    % plot eigenvalues to show fraction of variance captured
    figure
    plot(cumsum(eVals/sum(eVals)),'o','markeredgecolor',[1 0.5 0.1],...
        'markerfacecolor', [1 0.5 0.1],'markersize',8)
    %adjust font and font size
    set(gca, 'FontName', 'Helvetica', 'FontSize', 16);

    % plot the eigenworms
    figure
    for i = 1:numEigWorms
        subplot(ceil(numEigWorms/2),2,i)
        plot(eigenWorms(i,:), 'Color', [1 0.5 0.1], 'LineWidth', 2)
        xlim([0 size(eigenWorms, 2) + 1])
        title(num2str(i))
        %adjust font and font size
        set(gca, 'FontName', 'Helvetica', 'FontSize', 13);
    end
    
    % plot the covariance matrix
    figure
    imagesc(covarianceMat)
    set(gca,'YDir','normal')
    %adjust font and font size
    set(gca, 'FontName', 'Helvetica', 'FontSize', 16);
end