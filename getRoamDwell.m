function roamdwell = getRoamDwell( angspeed, speed, onfood, intlen )
%GETROAMDWELL.m This function annotates worm movement trajectories on food
%for roaming and dwelling
speed(~onfood) = NaN;
angspeed(~onfood) = NaN;
angspeed = abs(angspeed);

angspeed_intmean = movmean(angspeed,intlen)';
speed_intmean = movmean(speed,intlen)';
goodinds = ~isnan(angspeed_intmean) & ~isnan(speed_intmean);
angspeed_intmean = angspeed_intmean(goodinds);
speed_intmean = speed_intmean(goodinds);

end

