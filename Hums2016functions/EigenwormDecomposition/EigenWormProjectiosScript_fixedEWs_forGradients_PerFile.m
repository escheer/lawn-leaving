%This script projects anlge time series onto fixed input Eigenworms
%and saves data for each experiment into structure 'EigenWormTracks'
% consisting of undulation mode (termed 'CBW' here) and turning mode
% (termed 'RBW' here), full body posture (termed 'BW' here)
%
%(c) Ingrid Hums, ingrid.hums@gmail.com
%Created 2015
%
%code dependencies:
% EigenWormDecomp_fixedEWs



clearvars -except Eigenworms

files = dir('*_anglesV8*-corrected.mat');
   
    

for exp = 1: length(files)
    
    clear EigenWormTracks RBW CBW Tracks
    
    SaveName = [files(exp).name(1:end-31) 'V8_EigenWormTracks.mat'];
    
    load(files(exp).name);
    load(files(exp).name(1:end-32));
    

     for CW = 1:length(WormGeomTracks) 

            [CBW, RBW] = EigenWormDecomp_fixedEWs(WormGeomTracks(CW).WormAnglesNewRAD,Eigenworms);

             EigenWormTracks(CW).CBW = CBW';
             EigenWormTracks(CW).RBW = RBW';

             EigenWormTracks(CW).BW = WormGeomTracks(CW).WormAnglesNewRAD;
                
     end

        save(SaveName,'EigenWormTracks');
   
end

