function [projectionsEW12,projectionsEW3end,coefs,scoresAcc,Variances] = EigenWormDecomp(allangles)
%This function calculates the eigenworms (EWs)
%and reconstructs undulation mode (EW1 + EW2) and 
%turning mode (SUM of EW3-AngNum) of input angle time series
%by projecting angle data onto the eigenworms.
%Rows of X correspond to observations, columns to variables.
%
%(c) Ingrid Hums, ingrid.hums@gmail.com
%Created 2015
%
%Needs "princomp" from MatLab statistics toolbox.
%
%code dependencies:
% fixnanC



NumFrames = size(allangles,1);
AngNum = size(allangles,2);

projectionsEW12 = NaN(AngNum,NumFrames,'single')';
projectionsEW3end = NaN(AngNum,NumFrames,'single')';


msgid = 'stats:princomp:colRankDefX';
 warning('off', msgid);
 
         
        [coefs,scoresAcc,Variances,~] = princomp(fixnanC(allangles));
        % PRINCOMP assumes: Rows of X correspond to observations, columns to variables
        % COEFF is a P-by-P matrix, each column containing coefficients
        % for one principal component.  The columns are in order of decreasing
        % component variance. Rows of SCORE correspond to observations, columns to components.
        % coefs = PCs or eigenvectors (columns)
        % scoresAcc = projection amplitude (columns), representation of X in the principal component space
        % variances = eigenvalues of covariances matrix

        
      % Undulation Mode is made up by first 2 EW (independent of total number of EW) 
        
        finaldata = [coefs(:,1) coefs(:,2)]' * (fixnanC(allangles) - repmat(nanmean(allangles),NumFrames,1))';
        projectionsEW12(:,:) = finaldata' * [coefs(:,1) coefs(:,2)]' + repmat(nanmean(allangles),NumFrames,1);
        projectionsEW12 = projectionsEW12';
        projectionsEW12(isnan(allangles)') = NaN;


      % Turning Mode is made up from all remaining EW (or can possibly be restricted
      % to 3-4 or 3-9 or 3-whatever, as the lower EW have very little eigenvalues)

       featurevector = [];
       for n = 3: AngNum   
           featurevector = [featurevector coefs(:,n)];
       end
       
        finaldata = featurevector' * (fixnanC(allangles) - repmat(nanmean(allangles),NumFrames,1))';
        projectionsEW3end(:,:) = finaldata' * featurevector' + repmat(nanmean(allangles),NumFrames,1);
        projectionsEW3end = projectionsEW3end';
        projectionsEW3end(isnan(allangles)') = NaN;
        
    end
    


 