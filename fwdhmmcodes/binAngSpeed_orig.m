function binnedAngSpeed = binAngSpeed(trackData, binSize) %binSize is in frames
	for (i=1:length(trackData)) 
        angspeedData = trackData(i).AngSpeed;
        numbBins = (length(angspeedData))/binSize;
        for (j=1:numbBins)
            startInd = (j*binSize) - (binSize -1);
            endInd = (j*binSize);
            currentData = angspeedData(startInd:endInd);
            binnedAngSpeed(i).AngSpeed(j) = mean(abs(currentData));
        end
    end
    
end
