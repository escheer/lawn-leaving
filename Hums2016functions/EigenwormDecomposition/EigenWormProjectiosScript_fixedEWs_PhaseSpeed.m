%This script caclulates phasespeed = undulation frequency:
%It projects angle time series onto fixed input eigenworms EW1 & EW2
%to retrieve respective projection amplitudes and calculates the 
%phase speed from that for each worm track.
%
%(c) Ingrid Hums, ingrid.hums@gmail.com & Manuel Zimmer, manuel.zimmer@imp.ac.at
%Created 2015
%
% requires MatLab function smooth from Curve Fitting Toolbox)
%
%code dependencies:
% AccRevDatsV2b
% EigenWormDecomp_fixedEWs_Projections12




Projections = [];

SaveMe = 0;
SaveName = 'test';

framerate = 10; % fps
smoothfactor = 5;


if ~exist('WormGeomTracks','var')
    [WormGeomTracks, files, DatasetPointer] = AccRevDatsV2b('*_anglesV8-corrected.mat'); 
end

try
    WormGeomTracks = rmfield(WormGeomTracks,{'WormPointsNew','WormSegmentVectorsNew','WormHeadPosition','Status','SignSwitchedFrames'});
end


NumTracks = length(WormGeomTracks);

for CW = 1:NumTracks  
    
    [PA] = EigenWormDecomp_fixedEWs_Projections12(WormGeomTracks(CW).WormAnglesNewRAD,Eigenworms);
    Projections(CW).EW12 = PA;

    NaNpointer = isnan(sum(WormGeomTracks(CW).WormAnglesNewRAD,2));

    PC1 = smooth(Projections(CW).EW12(:,1),smoothfactor,'lowess');
    PC2 = smooth(Projections(CW).EW12(:,2),smoothfactor,'lowess');

    PC1(NaNpointer) = NaN;
    PC2(NaNpointer) = NaN;
    
    Projections(CW).EW12smooth(:,1) = PC1;
    Projections(CW).EW12smooth(:,2) = PC2;
  
    PhaseAngle = NaN(size(WormGeomTracks(CW).(AngAttr),1),1);
    
    for CF = 2:size(WormGeomTracks(CW).(AngAttr),1)  % get the phase angle by calculating the cross product of adjacent vectors in PC1 PC2 space

        vec1 = [0;PC1(CF);PC2(CF)];
        vec2 = [0;PC1(CF-1);PC2(CF-1)];

        CP = cross(vec2/norm(vec2),vec1/norm(vec1)); 
        PhaseVec = CP(1); 
        PhaseAngle(CF) = asin(PhaseVec);

    end
    
    Projections(CW).PhaseSpeed = PhaseAngle * framerate / (2*pi); % cycles per second
        
end
 

if SaveMe
    
    save([SaveName '-Projections'],'Projections','SaveName');

end




