function binnedSpeed = binSpeed(trackData, binSize) % binSize is in Frames here
	for (i=1:length(trackData)) 
        speedData = trackData(i).Speed;
        numbBins = (length(speedData))/binSize;
        for (j=1:numbBins)
            startInd = (j*binSize) - (binSize -1);
            endInd = (j*binSize);
            currentData = speedData(startInd:endInd);
            binnedSpeed(i).Speed(j) = mean(currentData);
        end
    end
end
