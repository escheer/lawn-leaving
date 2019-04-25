function omat=fixnanR(mat)
%omat=fixnan(mat)
%remove NaN values from matrix by selecting value above NaN areas within ROW
%(assumes ROWS are timeseries) !!!!
%(c) Saul Kato, saul@kato.com
%Created 121101
%
%updated 1301028 to check first row for NaNs and make them zeros and do NaN blocks
% Revised by Ingrid Hums to fix rows


%fill in all nans
for i=1:size(mat,1)
    if isnan(mat(i,1)); mat(i,1)=0; end
    for t=1:size(mat,2)
        if isnan(mat(i,t)); mat(i,t)=mat(i,t-1); end
    end
end
omat=mat;
end