% needs Tracks_als files as input.
% plots all tracks from 1 behavioral experiment at a time and lets the user define roi by clicking
%
%(c) Julia Riedl, julia.riedl@imp.ac.at
%Created 2016


files=dir('*als.mat');
roi={NaN(2,2)};


for F=1:length(files)
    
    load(files(F).name)
    figure('Visible','on')
    
    for i= 1:2:length(Tracks)
        %plot the track and starting point:
        hold on
        
        if length(Tracks(i).SmoothX)>150
            if F<12
                plot(Tracks(1,i).SmoothY,Tracks(1,i).SmoothX,'k')
            else
                plot(Tracks(1,i).SmoothX,Tracks(1,i).SmoothY,'k')
            end
            
        end
        
    end
    title(F)
    axis image
    pause (1)
    roi{F}(1,:)=ginput(1);
    roi{F}(2,:)=ginput(1);
    
    close
end

if ~isempty(files)
    save roinew roi
end




