function statemap = getStateAuto(trackData,binSize,ratio) % binSize in Frames here
binnedSpeed = binSpeed(trackData,binSize);
binnedAngSpeed = binAngSpeed(trackData,binSize);
statemap = getState(binnedSpeed,binnedAngSpeed,ratio); % 450 determined empirically
end
