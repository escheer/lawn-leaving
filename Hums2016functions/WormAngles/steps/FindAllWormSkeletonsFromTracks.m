%This script calculates skeletons (output: WormTraces) from binary worm images. 
%It requires the Parallel WormTracker format with the
%additional modifciation that each Tracks structure entry contains a
%binary worm image centered around the worm's centroid coordinates: Tracks.WormImages{}
%The script runs through all tracks files in the current directory.
%
%(c) Manuel Zimmer, manuel.zimmer@imp.ac.at
%Created 2103
%Revised 2015 by Ingrid Hums, ingrid.hums@gmail.com
%
% requires MatLab Image Processing Toolbox
%
%code dependencies:
% Skel2Line


close all;
clear;

FilenameString = '*tracks_als.mat';
files = dir(FilenameString);
[NumFiles,~] = size(files);



for CurrentFile = 1:NumFiles
    
    tic;
    
    close all;
    
    FileToAnalyze = files(CurrentFile).name
    
    load(FileToAnalyze);

    WormGeomTracks = struct('WormTraces', repmat({{}},1,length(Tracks)) );
    
    [~, NumTracks] = size(Tracks);
        
    [~, FileBaseName, ~] = fileparts(FileToAnalyze)
   
    length(Tracks)
    
    for  Index = 1:length(Tracks)
        
        CurrentTrack = (Index);
        
        disp([FileToAnalyze ':' 'current track: ' num2str(CurrentTrack)]);
                
        WormTraces = {};
        
        WormSegmentVectors = [];
        
        NumFrames = Tracks(CurrentTrack).NumFrames;
        
        for CurrentFrame=1:NumFrames
            
            close all
            
            CurrentFrameAlsStatus=1;
                        
            fr=Tracks(CurrentTrack).WormImages{CurrentFrame};
            
            [WormImageSizeY,  WormImageSizeX] = size(fr);
            
            fr2 = zeros(WormImageSizeY+10, WormImageSizeX+10);
            
            fr2(6:end-5,6:end-5) = fr;
            
            bwim = ~im2bw(fr2,0.1); % figure;imshow(bwim); %0.1
            
            bwimcomp = imcomplement(bwim); % switch black and white
            img1 = imdilate(bwimcomp, strel('disk',1)); %figure;imshow(img); % dilate parts of worm not fully filled (e.g. head sometimes)
             %img2 = bwmorph(img1,'erode',.5); 
            img2 = imerode(img1, strel('disk', 1));
            img3 = imcomplement(img2);   % figure; imshow(img3);  % switch b-w back
            
            edgewm = edge(img3,'canny',0.5,1.8); % figure; imshow(edgewm); % 0.05,3 original;  NEW: use 0.05,1 
            
            bridgedworm = bwmorph(edgewm,'bridge'); % figure; imshow(bridgedworm);
            
            filledworm = imfill(bridgedworm,'holes'); % figure; imshow(filledworm);            
            
            erodedworm = bwmorph(filledworm,'erode',1);  % figure; imshow(erodedworm); % 1.5 original; NEW: 1
                      
            wormskeleton = bwmorph(erodedworm,'skel',inf); 
            
            [wormskeleton, RecursionDepth] = Skel2Line(wormskeleton,0,10);
%              figure('Position',[858-210 824 204 128]); imshow(wormskeleton);
            
            Endpoints = bwmorph(wormskeleton, 'endpoints');
            
            [EndpontIndcsX, EndpontIndcsY] = find(Endpoints==1);
            
            if EndpontIndcsX
                
                wrm = bwtraceboundary(wormskeleton, [EndpontIndcsX(1) EndpontIndcsY(1)], 'NE');

                %partition into subsegments
                if (wrm)
                    
                    wrmtr=wrm(1:(end+1)/2,:);
                    smoothwormX = wrmtr(:,2)';
                    smoothwormY = wrmtr(:,1)';
                    
                    if numel(smoothwormY) > 45; %commands below will cause a crash when skeleton shrinks too much
                        WormTraces{CurrentFrame} = [smoothwormX; smoothwormY]';
                    else
                        CurrentFrameAlsStatus = 0;
                    end
                else
                    CurrentFrameAlsStatus = 0; 
                end
            else
                CurrentFrameAlsStatus = 0; 
            end

            if     CurrentFrameAlsStatus == 0;
                WormTraces{CurrentFrame} = {};  
            end
        end
  
        WormGeomTracks(CurrentTrack).WormTraces = WormTraces;
   
    end % end Track loop
    
    toc;
    
    save([FileBaseName '_anglesV8'],'WormGeomTracks');

end

