function [Ampl] = CaclulateAmplitude(EigenWormTracks,Attr,TmpTrcks)
%amplitde is the sum of the absolutes of all segment angles per frame
%output is a 2D matrix (tracks x frames)
%
%(c) Ingrid Hums, ingrid.hums@gmail.com 
%Created 2015


NumTracks = length(EigenWormTracks);
NumFrames = max([TmpTrcks.Frames]);

Ampl = NaN(NumTracks,NumFrames,'single');


for CW = 1:length(EigenWormTracks)
    
    
    time = TmpTrcks(CW).Frames;
         
    UM = EigenWormTracks(CW).(Attr);
   
    Ampl(CW,time) = sum(abs(UM),2);
  
end


end  