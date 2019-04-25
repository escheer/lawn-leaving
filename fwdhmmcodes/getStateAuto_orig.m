
function statemap = getStateAuto(trackData,binSize) % binSize in Frames here
    binnedSpeed = binSpeed(trackData,binSize);
    binnedAngSpeed = binAngSpeed(trackData,binSize)
    statemap = getState(binnedSpeed,binnedAngSpeed,450) % 450 determined empirically
end
