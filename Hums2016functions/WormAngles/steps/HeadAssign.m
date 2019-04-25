%% assigns head position and writes it into skel matfile:
%
%(c) Julia Riedl, julia.riedl@imp.ac.at
%Created 2015
%
%code dependencies:
% FindHead


clear

flipvalue = 75;


% global rois
% cd('mat')

id = strfind(cd, '\');
% %find and load als and skeleton files
home = cd;


skelfilescorr = dir('*anglesV8*-corrected.mat');

file = dir('*tracks_als.mat');

if size(skelfilescorr,1) == size(file,1)


for F = 1:length(skelfilescorr)
    
    Last=0;
    First=0;
    Xlast=[];
    Ylast=[];
    
%     file=dir([skelfiles(F).name(1:end-13) '.mat']);
    
    load (file(F).name);
    
    
    load (skelfilescorr(F).name);
    
    display(['...head assignment:'  skelfilescorr(F).name 'n\' file(F).name])
    fc=2;  
    
    for track= 1: length(WormGeomTracksC) 
        
    Last=0;
    First=0;
    Xlast=[];
    Ylast=[];
        
        % get skeleton:
        Skeleton = WormGeomTracksC(1,track).WormPointsNew;
        CurrentTrack = Tracks(1,track);
        
        headX = NaN(1,6);
        headY = NaN(1,6);
        flips = 1;
        cc = 1;
        breakflag = 0;
        CurrentFrame = 1;

 
while CurrentFrame<length(CurrentTrack.Frames) % for breaking the loop if flip episodes grow too long
    
    flips=1;
    
  
    for CurrentFrame = cc:length(CurrentTrack.Frames)
        
        if mod(CurrentFrame,135)==0
            p = 1;
        end
        
        prevLast=Last;
        prevFirst=First;
        


        
        [headX(CurrentFrame), headY(CurrentFrame),First,Last,CurrentTrack,flips,cc,] = FindHead(CurrentTrack,CurrentFrame,Skeleton,prevLast,prevFirst,Xlast,Ylast,headX,flips,cc,breakflag,flipvalue);
        
        Xlast=headX(CurrentFrame);
        Ylast=headY(CurrentFrame);
        breakflag=0;
       
        % switch back if alternative head position was used too often (>75 flips)
        if flips > flipvalue
            fc=fc+1;
            frameflag(fc)=track;
            prevLast=First;
            prevFirst=Last;
            breakflag=1;
            display(['bf:' num2str(track)]);
            if frameflag(fc)==frameflag(fc-2)
                breakflag=0;
                flips=1;
            else
            break % break the loop and restart at cc with switched head position
            end
        end
        
    end
end
        
        WormGeomTracksC(1,track).WormHeadPosition = horzcat(headX', headY');
        
                
    end % end track loop
    
    
    WormGeomTracksC(1).Status = 'headassigned';
    WormGeomTracks = WormGeomTracksC;
    
    clear WormGeomTracksC
    
    save(skelfilescorr(F).name, 'WormGeomTracks');
  
end %end file loop


else
      display('...files missing')
      
end
