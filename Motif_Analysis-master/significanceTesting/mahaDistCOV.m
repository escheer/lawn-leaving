function distance = mahaDistCOV(cov1, cov2, n1, n2, xbar1, xbar2,...
    equalVariance)

% MAHADISTSHRINK  Calculate the Mahalanobis distance between two data
% populations given the covariance matrix of each data set.
%
% Input
%   cov1 and cov2   - covariance matrices of the data to be compared.
%   n1 and n2       - the number of samples in each data set that went into
%                     computing the two covariance matrices.
%   xbar1 and xbar2 - row vectors with the means for each feature of the
%                     data sets to compare.
%   equalVariance   - Logical.  Are the variances of the two data sets
%                     equal?
%
% Output
%   distance        - the two-sample Mahalanobis distance between the data
%                     sets
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


% get the means and mean differences of each feature for each data set
xdiff = xbar1 - xbar2;

% calculate the covariance matrix for equal or unequal variances
if (equalVariance)
    S = ((n1 - 1) * cov1 + (n2 - 1) * cov2)/(n1 + n2 - 2);
else
    S = cov1/n1 + cov2/n2;
end

% check the condition of the matrix S.  Don't bother calculating distance
% if matrix is singular.  Return NaN instead.
if rcond(S) > eps
    % calculate the Mahalanobis distance
    distance = sqrt(xdiff / S * xdiff');
else
    distance = NaN;
end