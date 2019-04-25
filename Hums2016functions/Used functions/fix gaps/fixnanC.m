function omat=fixnanC(mat)
%omat=fixnan(mat)
%remove NaN values from matrix by selecting value above NaN areas within column
%(assumes COLUMNS are timeseries) !!!!
%(c) Saul Kato, saul@kato.com
%Created 121101
%
%updated 1301028 to check first row for NaNs and make them zeros and do NaN blocks
% Revised by Ingrid Hums to fix columns


%fill in all nans
for i=1:size(mat,1)
    for t=1:size(mat,2)       
        if isnan(mat(1,t)); mat(1,t)=0; end          
        if isnan(mat(i,t)); mat(i,t)=mat(i-1,t); end
    end
end
omat=mat;
end