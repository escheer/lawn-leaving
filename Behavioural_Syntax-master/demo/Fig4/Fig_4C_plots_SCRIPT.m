% script to find the n-grams that appear as hits the most frequently across
% all the strains/conditions

n = 1;

load(['Fig_4C_WORKSPACE_' num2str(n) '-gram.mat'])

% get all of the hits and the corresponding comparison stats
allHits = vertcat(informativeBehaviours{:});
comparisonStatMat = vertcat(comparisonStats{:});


% count the hits
[uniqueHits, counts] = countUniqueRows(allHits);

% get the frequency rank of each of the unique hits in the different
% strains
rankCell = cell(size(uniqueHits, 1), 1);
for ii = 1:size(uniqueHits, 1)
    % find all examples of the current hit
    matchInds = ismember(allHits, uniqueHits(ii, :), 'rows');
    
    % get the rank of the match
    rankCell{ii} = comparisonStatMat(matchInds, 3);
end

% calculate the histogram of hits
edges = exp(1:log(5e3)/60:log(5e3));
rankHist = histc(vertcat(rankCell{:}), edges);
stairs(edges, rankHist, 'b')
set(gca, 'XScale', 'log')
set(gca,'ticklength',2*get(gca,'ticklength'))
saveas(gca, 'rankDistribution.eps')
    
% figure
% edges = 1:5e3/500:5e3;
% rankHist = histc(vertcat(rankCell{:}), edges);
% plot(edges, rankHist)
% set(gca, 'XScale', 'log')


% also make box plots for each of the strains for the different hits
for ii = 1:size(uniqueHits, 1)
    % get the index of the current hit in the total set of unique n-grams
    matchInd = ismember(uniqueNGrams, uniqueHits(ii, :), 'rows');
    
    % use anova1 function to plot box plot of frequencies in all strains
    [~, table, stats] = anova1(freqMat(:, matchInd), wormNames);
    
    
    % recalculate p-values for this n-gram to add bars to the plot
    sigDiffMap = zeros(numel(uniqueNames));
    offset = 5;
    for jj = 1:numel(uniqueNames)-1
        % get name matches of jj'th strain
        nameMatches1 = find(strcmp(wormNames, uniqueNames{jj}));
        freq1 = freqMat(nameMatches1, matchInd);
        
        for kk = jj+1:numel(uniqueNames)
            % get the submatrix from count mat for the current pair
            nameMatches2 = find(strcmp(wormNames, uniqueNames{kk}));
            freq2 = freqMat(nameMatches2, matchInd);
            
            % calculate the p-value
            p = ranksum(freq1, freq2);
            
            % for hits, add a bar to the box plot to indicate significant
            % difference
            if p <= pCrit
                % add the significantly different points to the sigDiffMap
                sigDiffMap(jj, kk) = 1;
                sigDiffMap(kk, jj) = 1;
                
                % add bars and p-values
                ylim([-5e-4, (offset+0.5) * stats.s])
                line([jj, kk], [offset * stats.s, offset * stats.s], ...
                    'LineWidth', 2, 'Color', [0.1 0.3 0.9])
                text((kk+jj)/2-0.25, (offset + 0.25) * stats.s, num2str(p), 'FontSize', 14)
                offset = offset + 0.5;
            end
        end
    end
    
    
    title(['Number of total hits: ' num2str(counts(ii)) '.'], ...
        'FontSize', 14)
    
    % save the box plot
    saveas(gcf, ['./boxplots-all-' num2str(n) ...
        '-gram/boxplot_' num2str(ii) '.pdf'])
    close all
    
    % plot and save the map of significant differences
    imagesc(sigDiffMap)
    colormap([1,1,1; 0,0,0])
    set(gca,'XTick',0.5:numel(uniqueNames)+0.5)
    set(gca,'YTick',0.5:numel(uniqueNames)+0.5)
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'GridLineStyle', '-')
    grid on
    
    saveas(gcf, ['./boxplots-all-' num2str(n) ...
        '-gram/sigDiffMap_' num2str(ii) '.pdf'])
    close all
    
    
    % also plot the behaviour
    figure
    plotSequence(postures, uniqueHits(ii, :), 'r', 'b')
    text(0, 0.6, num2str(uniqueHits(ii, :)), 'FontSize', 16)
    
    % save the behaviour plot
    saveas(gcf, ['./boxplots-all-' num2str(n) ...
        '-gram/posturePlot_' num2str(ii) '.pdf'])
    close all
end




