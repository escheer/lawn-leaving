function [TrcksAcc, files, DatasetPointer] = AccRevDatsV2(flnmstr)
%loads tracks data from wormanalyzer 'Analyze all tracks' into one
%structure 'TrcksAcc' and returns a cell array that contains all filenames
%(c) Manuel Zimmer, manuel.zimmer@imp.ac.at



flnms=dir(flnmstr); %create structure from filenames

files={flnms.name}';

[max, ~] = size(flnms);

TrcksAcc =[]; %structure that accumulates all tracks data

DatasetPointer(1,1) = 1;
DatasetPointer(1,2) = 0;

files{1,2} = 1;
files{2,3} = 0;

for i=1:max;
    
    i
   
    file = load (flnms(i).name);
    
    [~,NumTracks] = size(file.Tracks);
    
    if i>1
        
       DatasetPointer(i,1) = DatasetPointer(i-1,2)+1; 
       DatasetPointer(i,2) =  DatasetPointer(i-1,2) + NumTracks;
       
       files{i,2} = files{i-1,3}+1; 
       files{i,3} =  files{i-1,3} + NumTracks;
         
    else
        DatasetPointer(i,2) = NumTracks;
                files{i,3} = NumTracks;
        
    end
    
    TrcksAcc = [TrcksAcc file.Tracks];
    
end;