function statemap = getState(SpeedBins,AngSpeedBins,ratio) % converts sequence of Speed/AngSpeed into a binary of 1's (dwell bins) and 2's (roam bins)
    
    for(j=1:length(SpeedBins))
        nBins = length(SpeedBins(j).Speed);
        statemap(j).state = [];
        for (i=1:nBins)
            actualRatio = ((AngSpeedBins(j).AngSpeed(i))/(SpeedBins(j).Speed(i)));
            if((isnan(actualRatio)) == 1)
                if(i==1)
                    statemap(j).state(i) = NaN;
                else
                statemap(j).state(i) = statemap(j).state(i-1);
                end
            else
            if (actualRatio >= ratio)
                statemap(j).state(i) = 1;
            else
                statemap(j).state(i) = 2;
            end
            end
        end
    end
    
end
