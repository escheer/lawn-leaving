function bindata = meanbindata( timeseries,binwidth,fn )
%BINDATA.m This function takes a timeseries and applies a mean to each FULL
%bin (partially filled bins omitted). nanmean is used so nans are ok.
[bins,edges] = discretize(1:length(timeseries),1:binwidth:length(timeseries));

bindata = zeros(1,length(edges)-1);
for i = 1:length(edges)-1
    bindata(i) = nanmean(timeseries(bins==i));
end

end

