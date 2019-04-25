function distance = mahaDistShrink(data1, data2, equalVariance)

% MAHADISTSHRINK  Calculate the Mahalanobis distance between two data
% populations using the James-Stein type shrinkage estimator of the
% covariance matrix to deal with having a larger number of features than
% number of samples.
%
% Input
%   data1 and data2 - two data sets of feature vectors to be compared. The
%                     samples are arranged along rows with features along
%                     the columns
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


% get the number of samples in each data set
n1 = size(data1, 1);
n2 = size(data2, 1);

% get the shrinkage estimate of the covariance matrix of each data set
[S1, ~, ~] = covshrinkKPM(data1, 0);
[S2, ~, ~] = covshrinkKPM(data2, 0);

% get the means and mean differences of each feature for each data set
xbar1 = mean(data1);
xbar2 = mean(data2);
xdiff = xbar1 - xbar2;

% calculate the covariance matrix for equal or unequal variances
if (equalVariance)
    S = ((n1 - 1) * S1 + (n2 - 1) * S2)/(n1 + n2 - 2);
else
    S = S1/n1 + S2/n2;
end

% check the condition of the matrix S.  Don't bother calculating distance
% if matrix is singular.  Return NaN instead.
if rcond(S) > eps
    % calculate the Mahalanobis distance
    distance = sqrt(xdiff / S * xdiff');
else
    distance = NaN;
end