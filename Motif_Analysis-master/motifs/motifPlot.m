function motifPlot(motifDictionary, numDimensions, maxLength, fps, ...
    pauseLength, waitBetweenPlots)

% MOTIFPLOT  Make a series of plots (in a single plot window) showing the
% different motifs in a motif dictionary
% 
% Input
%   motifDictionary  - A cell array containing the parent file names, the
%                      actual motifs, and their start locations at the
%                      original time series resolution (i.e. no longer
%                      down sampled)
%   numDimensions    - The number of eigen worms to consider in motif
%                      finding
%   maxLength        - The longest motif in the dictionary
%   fps              - The frame rate to convert frames to seconds for plot
%   pauseLength      - How long should the pause be between motif plots?
%   waitBetweenPlots - Should MATLAB wait for the user to press a button
%                      between motif plots to allow for inspection?
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

h = 1;
for i = 1:size(motifDictionary, 1)
    % should each plot be displayed until a keypress?
    if waitBetweenPlots
        waitforbuttonpress
    end
    % check that the current entry is not empty
    if ~isempty(motifDictionary{i, 1})
        
        hold off
        % plot each dimension separately
        for j = 1:numDimensions
            % add a zero line
            plot((-5:length(motifDictionary{i,3}(j, :))) / fps,...
                10*j, 'Color', [0.6 0.6 0.6]);
            hold on
            
            % plot the amplitudes
            plot((1:length(motifDictionary{i,3}(j, :))) / fps,...
                motifDictionary{i,3}(j, :) + 10*j, 'Color', [1 0.6 0.3],...
                'LineWidth', 1.1)
            plot((1:length(motifDictionary{i,5}(j, :))) / fps,...
                motifDictionary{i,5}(j, :) + 10*j, 'Color', [0.3 0.6 1],...
                'LineWidth', 1.1)
        end
        % set plot dimensions, get frame
        xlim([0 maxLength] / fps)
        ylim([0 10*numDimensions + 10])
        
        % add text annotation indicating position in motifIndex
        text(maxLength * 0.9, (10*numDimensions + 10) * 0.9, num2str(i))
        
        getframe(h);
        
        % wait pauseLength seconds to allow for viewing
        pause(pauseLength)
    end
end
