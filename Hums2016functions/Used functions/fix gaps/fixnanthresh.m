function [xnew] = fixnanthresh(x,thr)
% fix NaNs stretches by replacing with previous value in series
% only NaN stretches <= threshold length (thr)
% works only for a vector;

nidx = isnan(x(:))';
nstart = strfind([0, nidx], [0 1]);
nlen = strfind([nidx, 0], [1 0]) - nstart + 1;

%find indices of NaN elements (won't be interpolated)
nstart = nstart(nlen > thr);
nend = nstart + nlen(nlen > thr) - 1;
idx = cell2mat(arrayfun(@colon, nstart, nend, 'UniformOutput', false));

%interpolate and replace original NaN elements with NaN values
 xnew = x;
 
 if length(xnew(~nidx)) > 1 
    if isnan(xnew(1)); xnew(1)=0; end
    for t=1:length(x)
        if isnan(xnew(t)); xnew(t)=xnew(t-1); end
    end
    xnew(idx) = NaN;  
 end


end

