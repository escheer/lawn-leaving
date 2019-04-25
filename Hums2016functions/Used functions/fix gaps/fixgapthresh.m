function [xnew] = fixgapthresh(x,thr,inter)
% interpolate over NaNs, only NaN stretches <= threshold length (thr)
% with interpolation method specified by input 'inter'
% works only for a vector (Row or Column, but output is ROW)

nidx = isnan(x(:))';
nstart = strfind([0, nidx], [0 1]);
nlen = strfind([nidx, 0], [1 0]) - nstart + 1;

%find indices of NaN elements (won't be interpolated)
nstart = nstart(nlen > thr);
nend = nstart + nlen(nlen > thr) - 1;
idx = cell2mat(arrayfun(@colon, nstart, nend, 'UniformOutput', false));

%interpolate and replace original NaN elements with NaN values
if length(x(~nidx)) > 1
    xnew = interp1(find(~nidx), x(~nidx), 1:numel(x), inter);
    xnew(idx) = NaN;
else
    xnew = x;
end

end

