%(1) reorders worm skeleton based on correct head assignment:
%flips skeleton order when first point is not equal to head coordinates
%(2) corrects left-right flips of the data ()
%
%(c) Julia Riedl, julia.riedl@imp.ac.at & Ingrid Hums, ingrid.hums@gmail.com
%Created 2015
% 
%code dependencies:
% fixnanthresh
% fixgapthresh



threshes = [1.5 3];
% Left-right flips are checked in two rounds with increasing threshold.
% You should change these threshes if there are still many LEFT-RIGHT flips!
% Checks the sum of the difference of segment angles between adjacent frames
% l-r flips are deteced when this sum is higher than the threshold
% i.e. posture cannot change that much between two frames (0.1 second)!
% second round switches back falsely corrected flips



angle_files=dir('*anglesV8*-corrected.mat');

for F = 1:length(angle_files) %load data
    
    tic
    load (angle_files(F).name);

if ~strcmp(WormGeomTracks(1).Status,'ordered')
    
    display(['...running on: ' angle_files(F).name])
    
   
    numIP = size(WormGeomTracks(1).WormAnglesNewRAD,2)+2;
    
    for T = 1:length(WormGeomTracks)
        
        headXY = WormGeomTracks(T).WormHeadPosition';
        
        WormPointsOrdered = {NaN};
        WormAngles = NaN(length(headXY),numIP-2);
        WormVectors = NaN(length(headXY),2,numIP-1);
        
        
        for fr = 1:length(WormGeomTracks(T).WormPointsNew)
                        
               
                if ~isempty(WormGeomTracks(T).WormPointsNew{fr})
                    
                    skeleton = WormGeomTracks(T).WormPointsNew{fr};
                    angles   = WormGeomTracks(T).WormAnglesNewRAD(fr,:);
                    vectors  = squeeze(WormGeomTracks(T).WormSegmentVectorsNew(fr,:,:));
                                        
                    if ~isempty(skeleton) && ~all(skeleton(1,:)' == headXY(:,fr))                      
                        
                        WormPointsOrdered{fr} = flipud(skeleton);
                        WormAngles(fr,:) = fliplr(angles);
                        WormVectors(fr,:,:) = -fliplr(vectors);
                        
                    else
                        
                        WormPointsOrdered{fr} = skeleton;
                        WormAngles(fr,:) = angles;
                        WormVectors(fr,:,:) = vectors;
                        
                    end 
                    
                else
                    WormPointsOrdered{fr} = {};
                end
                
        end
        
        % jumps: probably when worm is not heading E or W but sharply N or S
        % angle is sign flipped
        % do this 2x, in 2nd round get some more or wrongly flipped with
        % higher threshold
        
       switches = []; % note done which frames were switched

        for h = 1:length(threshes)
            
            fixed = NaN(size(WormAngles));
            for a = 1:numIP-2
             fixed(:,a) = fixnanthresh(WormAngles(:,a),2); % small gaps are fixnaned so that flip afterwards can still be detected
             fixed(:,a) = fixgapthresh(fixed(:,a),10,'cubic'); % bigger gaps must be interpolated to not falsely detect flips
            end
            jumps = diff(fixed).^2; % square to better extract the jumps (and make all values positive!)
            jumps = [0; sum(jumps,2)];

            %%
%              figure('Position',[3 800 1915 109]); plot(jumps,'.'); title(num2str(h));
%              figure('Position',[3 519 1915 109]);imagesc(fixed'); title(num2str(h));
            %%

             jumps = find(jumps > threshes(h)); % IMPORTANT THRESHOLD

            if mod(length(jumps),2)~=0
                jumps(end+1)= size(WormGeomTracks(T).WormAnglesNewRAD,1)+1;
            end
            i = 1:length(jumps);

            evens = i(mod(i,2)==0); 
            odds = i(mod(i,2)~=0); 
            starts = jumps(odds);
            ends = jumps(evens)-1;
            
            changeframes = [];
            for i = 1:length(starts)
                changeframes = [changeframes,starts(i):ends(i)];
            end
            
            WormAngles(changeframes,:) = WormAngles(changeframes,:) * (-1);
              
            for i = 1:length(changeframes)  
                WormPointsOrdered{changeframes(i)} = fliplr(WormPointsOrdered{changeframes(i)});
             
                points = WormPointsOrdered{changeframes(i)};
                if ~isempty(points)
                    WormVectors(fr,1,:) = points(2:end,1)-points(1:end-1,1);
                    WormVectors(fr,2,:) = points(2:end,2)-points(1:end-1,2);
                else
                    WormVectors(fr,:,:) = NaN;
                end
            end
            %%
            
            sw = [starts, ends];
            switches = [switches; sw ; NaN NaN];  % note done which frames were switched
            
%           figure('Position',[3 519 1915 109]); imagesc(WormAngles');
        end
        
        
        WormGeomTracks(1,T).WormPointsNew = WormPointsOrdered;

        WormGeomTracks(1,T).WormAnglesNewRAD = WormAngles;
        
        WormGeomTracks(1,T).WormSegmentVectorsNew = WormVectors;
        
        WormGeomTracks(1,T).SignSwitchedFrames = switches;
        
        
    end % end tracks loop
    
    if isfield(WormGeomTracks,'GradientLandscape')
        try
            WormGeomTracks=rmfield(WormGeomTracks, {'ReversalState','GradientLandscape'});
        end
    end
    
    %save new:
    WormGeomTracks(1).Status = 'ordered';
    toc
    display('save...')
    save(angle_files(F).name,'WormGeomTracks');
    toc
    
end

end