%This script projects anlge time series onto fixed Eigenworms 
%(need to be loaded into in workspace)
%and saves data from all experiments into one big structure 'EigenWormTracks': 
%undulation mode + turning mode;
%saves undulation + turning amplitude, body amplitude (=sum of the two amplitudes)
%as 2D matrices (tracks x time);
%also saves all angles time series (= body posture) into one structure 'WormGeomTracksRAD;
%'SBinTrcksSpd': 2D matrix (tracks x time) of speed (fwd + rev) exlcuding
%data at the assay arena borders;
%'RevState': 2D matrix (tracks x time) 0=forward, 1=reverse, NaN=no data;
%'CopperRing': 2D logical matrix, 1=data close to border, 0=rest
%'DataToNaN': 2D logical matrix, combines RevState and CopperRing:
%0=forward, 1=reverse or data at border or no data
%TmpTrcks stores the frames of each track with respect to the recording 
%DatasetPointer stores start and end track of each experiment with respect
%to the 2D matrices containing all data
%
%(c) Ingrid Hums, ingrid.hums@gmail.com
%Created 2015
%
%code dependencies:
% AccRevDatsV2b
% EigenWormDecomp_fixedEWs
% CaclulateAmplitude
% AccRevDatsV2
% spdalsV5_MSv103_woPir
% ShallowDeepOTalsV4_MSv107_IH_woPir


%%%%%%%%%%%%%%%%%%%%
% decide whether to save and under which same & set your video framerate
SaveMe = 0;

% SaveName = 'test';

framerate = 10;
%%%%%%%%%%%%%%%%%%%%


if ~exist('WormGeomTracks','var')
    [WormGeomTracks, files, DatasetPointer] = AccRevDatsV2b('*_anglesV8*-corrected.mat'); 
end

try
    WormGeomTracks = rmfield(WormGeomTracks,{'WormPointsNew','WormSegmentVectorsNew','WormHeadPosition','Status','SignSwitchedFrames'});
end


NumTracks = length(WormGeomTracks);


if ~exist('TmpTrcks','var') 
    if ~exist('Tracks','var')
        [Tracks, files, ~] = AccRevDatsV2('*_tracks_als.mat');
    end
    
    for i = 1:NumTracks
        TmpTrcks(i).Frames = Tracks(i).Frames;
    end
end



for CW = 1:NumTracks   
        [UM, TM] = EigenWormDecomp_fixedEWs(WormGeomTracks(CW).WormAnglesNewRAD,Eigenworms);
         EigenWormTracks(CW).UndulationMode = UM';
         EigenWormTracks(CW).TurningMode = TM';
end



 UndulationAmplitude = CaclulateAmplitude(EigenWormTracks,'UndulationMode',TmpTrcks);
 TurningAmplitude = CaclulateAmplitude(EigenWormTracks,'TurningMode',TmpTrcks);
 BodyAmplitude = CaclulateAmplitude(WormGeomTracks,'WormAnglesNewRAD',TmpTrcks);

 
 
%%
if ~exist('DataToNaN','var')
    
    if ~exist('Tracks','var')
        [Tracks, files, ~] = AccRevDatsV2('*_tracks_als.mat');
    end
    
    sParam.spdalsRingLimit = 80; 
    [~,  SBinTrcksSpd] = spdalsV5_MSv103_woPir(sParam,Tracks,1,10); 
    % detects worm data close to CopperRing border
    
    [~, NumTracks] = size(Tracks);
    MaxNumFrames = max([Tracks.Frames]);

    DataToNaN = isnan(SBinTrcksSpd);
    CopperRing = isnan(SBinTrcksSpd);
 
    
     sParam.ROalsRingLimit = 80;
    [~, ~, ~, ~, BinTrcksLRstate, BinTrcksSRstate] = ShallowDeepOTalsV4_MSv107_IH_woPir(sParam, Tracks, 1, framerate, 1);
    % finds reversals

    DataToNaN( BinTrcksLRstate == 1 | BinTrcksSRstate == 1) = 1; 
    RevState = BinTrcksLRstate; RevState(BinTrcksSRstate == 1) = 1;

end

clear Tracks
 

if SaveMe
    save([SaveName '-EigenWormTracks'],'EigenWormTracks','SaveName','UndulationAmplitude','TurningAmplitude','BodyAmplitude');
    save([SaveName  '_WormGeomTracksRAD'], 'SaveName','WormGeomTracks','SBinTrcksSpd','DatasetPointer','files','TmpTrcks','DataToNaN','RevState','CopperRing');
end




