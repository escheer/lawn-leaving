function [t2 pVal] = hotellingT2(data1, data2, equalVariance, permutations)

% HOTELLINGT2 Implements a Hotelling's T^2 test that can handle large
% dimensional data (more dimensions than samples) based on a shrinkage 
% estimate of the covariance matrix.  The algorithm is described in 
% Chen-An Tsai and James J. Chen (2009) Multivariate analysis of 
% variance test for gene set analysis, Bioinformatics 25:897-903.
% 
% If there are sufficient samples, then the covariance matrix is computed
% directly.
%
% This code is based on their R implementation downloaded from:
% http://mail.cmu.edu.tw/~catsai/research.htm (Dec 2011)
%
% Input
%   data1 and data2 - Matrices of data with samples along rows and features
%                     along columns.
%   equalVariance   - Logical.  Are the variances of the two data sets
%                     equal?
%   permutations    - The number of permutations to perform in the
%                     permutation test.  Set to 0 to skip permutation test.
%                     In this case pVal will be returned as NaN.
%
% Output
%   t2              - The Hotelling T^2 statistic for the comparison of
%                     data1 and data2
%   pVal            - The p-value for the test calculated using a
%                     permutation test.
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

% calculate the T^2 statistic for the comparison of the two data sets
t2 = hott2(data1, data2, equalVariance);

% calculate p-value using a permutation test
% make a combined data matrix for easy permutation
combinedData = [data1; data2];

% initialise hit number (number of times the permuted data has a larger T^2
% statistic than the original data)
hits = 0;

for i = 1:permutations
    % get a random permutation of the indices
    permutationInds = randperm( (size(data1, 1) + size(data2, 1)) );
    
    % remake the permuted data sets
    data1Perm = combinedData( permutationInds(1:size(data1, 1)), :);
    data2Perm = combinedData( permutationInds(size(data1, 1) + 1:end), :);
    
    t2Perm = hott2(data1Perm, data2Perm, equalVariance);
    
    if t2Perm >= t2
        hits = hits + 1;
    end
end

% calculate the p-value
if permutations == 0
    pVal = NaN;
else
    pVal = hits/permutations;
end


%--------------------------------------------------------%
% The actual Hotelling T^2 function described/used above %
%--------------------------------------------------------%

function t2 = hott2(data1, data2, equalVariance)

% get the number of samples for each data set
nn = size(data1, 1);
nd = size(data2, 1);

% get the covariance matrix of each data set.  If the sample size is less
% than or equal to the number of dimensions, us a shrinkage estimate of the
% covariance matrix.
if nn <= size(data1, 2)
    [SN, ~, ~] = covshrinkKPM(data1, 0);
else
    SN = cov(data1);
end
if nd <= size(data2, 2);
    [SD, ~, ~] = covshrinkKPM(data2, 0);
else
    SD = cov(data2);
end

% get the means and mean differences of each feature for each data set
xbar1 = mean(data1);
xbar2 = mean(data2);
xdiff = xbar1 - xbar2;

% calculate the Hotelling's T^2 statistic appropriate for equal or unequal
% variances
if (equalVariance)
    S = ((nd - 1) * SD + (nn - 1) * SN)/(nd + nn - 2);
    t2 = ((nd * nn)/(nd + nn)) * (xdiff / S * xdiff');
else
    S = SN/nn + SD/nd;
    t2 = xdiff / S * xdiff';
end
