function SUMMARY_STRUCT = fixRoamingDwelling(SUMMARY_STRUCT,centSmthWindow,RD_struct,EXIT_STRUCT,POKE_STRUCT,stat_int,titlestr)
%fixRoamingDwelling This function edits a SUMMARY_STRUCT to recompute its
%ANGSPEED and re-computes Roaming and Dwelling by 2D and HMM methods.

trans = RD_struct.trans;
emis = RD_struct.emis;
cutoff = RD_struct.cutoff;
x_offset = RD_struct.x_offset;
binSize = RD_struct.binSize;

if centSmthWindow==1
    CENTROID_SMTH = SUMMARY_STRUCT.CENTROID;
else
    CENTROID_SMTH = [movmean(SUMMARY_STRUCT.CENTROID(:,1),centSmthWindow,'omitnan') movmean(SUMMARY_STRUCT.CENTROID(:,2),centSmthWindow,'omitnan')];
end
%  ANGSPEED = getAngularSpeed(CENTROID_SMTH(:,1), CENTROID_SMTH(:,2));
% [ANGSPEED, ~] = getAngularSpeed2(CENTROID_SMTH, angSpeedWindow);
% [ANGSPEED, ~] = getAngularSpeed3(CENTROID_SMTH, angSpeedWindow);
% [ANGSPEED, ~] = getAngularSpeed4(CENTROID_SMTH,angSpeedWindow);

ANGSPEED = getAngularSpeed_NavinMethod(CENTROID_SMTH);

[ROAMDWELL_2D, ROAMDWELL_HMM, ~, ~, ~] = getHMMStatesConcatenatedTracks(SUMMARY_STRUCT.SPEED,ANGSPEED,binSize,SUMMARY_STRUCT.CENTROIDINLAWN,cutoff,trans,emis,x_offset,false);

SUMMARY_STRUCT.CENTROID_SMTH;
SUMMARY_STRUCT.ANGSPEED = ANGSPEED;
SUMMARY_STRUCT.ROAMDWELL_2D = ROAMDWELL_2D;
SUMMARY_STRUCT.ROAMDWELL_HMM = ROAMDWELL_HMM;

FRAC_ROAMING_HMM = nansum(ROAMDWELL_HMM==2)/sum(~isnan(ROAMDWELL_HMM));
SUMMARY_STRUCT.FRAC_ROAMING_HMM = FRAC_ROAMING_HMM;

makeSummaryFig(SUMMARY_STRUCT,EXIT_STRUCT,POKE_STRUCT,stat_int,titlestr); %eventually need to fix this to give it an outpath

end

