%This script calculates Eigenworms, derived from all input anlge time series. 
%Therefore it concatenates all angle time series in order to perform PCA on
%one large angle time series from all combined tracks.
%
%(c) Ingrid Hums, ingrid.hums@gmail.com
%Created 2015
%
%code dependencies:
% AccRevDatsV2b
% EigenWormDecomp
%
% Needs "princomp" from MatLab statistics toolbox.


%%%%%%%%%%%%%%%%%%%%%%
% decide whether to save and under which same 

SaveMe = 0; % 1=yes

SaveName = 'test';

%%%%%%%%%%%%%%%%%%%%%%%



if ~exist('WormGeomTracks','var')
    [WormGeomTracks, files, DatasetPointer] = AccRevDatsV2b('*_anglesV8*-corrected.mat'); 
    WormGeomTracks = rmfield(WormGeomTracks,{'WormPointsNew','WormSegmentVectorsNew','WormHeadPosition','Status','SignSwitchedFrames'});
end


NumTracks = length(WormGeomTracks);

allangles = [];


for CW = 1:NumTracks
    
     c = size(allangles,1);

     allangles = [allangles; WormGeomTracks(CW).WormAnglesNewRAD];
    
end
    
    
[~,~,coefs,scores,variances] = EigenWormDecomp(allangles);
    
scores = scores';
      
    
% plot the first 9 PCs

 for pc = 1:9

   subplot(3,3,pc);
   plot([1:length(coefs)], coefs(:,pc));
   title(['eigenworm ' num2str(pc)],'Fontsize',10);
   ylim([-.75 .75]);
   xlim([1 length(coefs)]);
   
 end

 ha = axes('Position', [0 0 1 1], 'Xlim', [0 1], 'Ylim', [0 1], 'Box', 'off', ...
     'Visible', 'off', 'Units', 'normalized', 'clipping' , 'off');
    text(0.5, 1,['Eigenworms from all ' SaveName], 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 9, 'FontWeight', 'bold');
    

    
    Eigenworms = coefs;
    
    Variances = variances;
    RelVariances = variances / sum(variances);
    RelVariances(:,2) = cumsum(RelVariances);
    

if SaveMe
    save([SaveName '-EigenWorms'],'SaveName','Eigenworms','Variances','RelVariances');
end




