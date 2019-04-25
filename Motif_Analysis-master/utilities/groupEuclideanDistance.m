function distance = groupEuclideanDistance(data1, data2)

% GROUPEUCLIDEANDISTANCE  Calculates the distance between the medians of
% two data populations.  The median is taken for each feature value and so
% does not necessarily represent an actual member of the population
% directly.
%
% Input
%   data1 and data2 - two data sets of feature vectors to be compared. The
%                     samples are arranged along rows with features along
%                     the columns
%
% Output
%   distance        - the Euclidean distance between the medians of each
%                     feature in each population
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


% calculate the median feature vector representing each population
medianData1 = median(data1, 1);
medianData2 = median(data2, 1);

% calculate the distance between the median vectors
distance = sqrt(sum((medianData1 - medianData2).^2));