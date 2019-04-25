function [projectionsEW12] = EigenWormDecomp_fixedEWs_Projections12(allangles,coefs)
%This function performs eigenworm decomposition based on 
%INPUT EIGENWORMS.
%It returns the corresponding projection amplitudes of the 
%input angle time series projected onto eingenworms EW1 and EW2.
%
%(c) Ingrid Hums, ingrid.hums@gmail.com
%Created 2015
%
%code dependencies:
% fixnanC


NumFrames = size(allangles,1);

projectionsEW12 = NaN(2,NumFrames,'single')'; % is defined here to store data of type 'single', will be transposed later
 

finaldata = [coefs(:,1) coefs(:,2)]' * (fixnanC(allangles) - repmat(nanmean(allangles),NumFrames,1))';
projectionsEW12(:,:) = finaldata';
        
    
end
    


 