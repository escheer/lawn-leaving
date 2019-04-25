% use on revinfo files
% plotes reversal probability vs. bearing
%
%(c) Julia Riedl, julia.riedl@imp.ac.at
%Created 2016
%
%code dependencies:
% dirname2
%
% third party code dependencies: errorb


clear
warning off;
CC=1;
CM= jet(30);
leg=NaN;

for delay=6
    
    %plot?
    plotting=1;
    files =dir('*revinfo*');
    
    %parameters:
    
    first=0;
    last=180;
    binsize=20;
    
    
    
    % for each experiment: put all run data into one vector:
    mB= NaN(55,40);
    mC= NaN(55,40);
    mT= NaN(55,40);
    revN= NaN(55,40);
    revNnorm=NaN(55,40);
    c=0;
    
    for batch=1:length(files)
             
        load(files(batch).name);
        disp(files(batch).name);
        
        for F=1:length(bearing)
            
            if ~isempty(bearing{F})
                
                c=c+1;
                bearingAll=[];
                SpeedAll=[];
                revAll=[];
                
                
                for i=1:length(bearing{F})
                    if ~isempty(dBearing{F}{i})
                        bv=bearing{F}{i};
                        sp=speed_ctr{F}{i};
                        rev=reversals{F}{i};
                        
                        bearingAll=cat(2,bearingAll,bv);
                        SpeedAll=cat(2,SpeedAll,sp);
                        revAll=cat(2,revAll,rev);
                    end
                end
                
                bv=round(bearingAll*10)/10;
                
                %%%%%bin:
                cc=1;
                if ~isempty(bv)
                    
                    for i=first :binsize:(last-binsize)
                        
                        bin_idx= find(bv>i & bv<=i+binsize);
                        
                        bin_idx(bin_idx>length(bv)-delay)=[];
                        
                        if ~isempty(bin_idx)
                            
                            revV=revAll(bin_idx+delay);
                            binB=bv(bin_idx);
                            binS=SpeedAll(bin_idx);
                            %remove reversals which happen at very low speed episodes:
                            bi=find(binS(revV)<0.02);
                            ri=(find(revV));
                            revV(ri(bi))=0;
                            
                            mB(c,cc)=nanmean(binB);
                            revN(c,cc)=nansum(revV);
                            revNnorm(c,cc)=nansum(revV)/length(bin_idx);
                        end
                        
                        cc=cc+1;
                        
                    end
                    
                end
            end
        end
        
    end % end files loop
    
    revNnorm=revNnorm(1:c,1:cc-1);
    revN=revN(1:c,1:cc-1);
    mB=mB(1:c,1:cc-1);
    h=NaN;
    
    %% plot:(1) reversal rate binned by bearing
    name=dirname2(cd);
    
    if plotting==1
        
        sem=nanstd(revNnorm*3,1)/sqrt(c);
        hold on
        h(CC)=plot(nanmean(revNnorm*3,1),'color',CM(delay,:));
        errorb(nanmean(revNnorm*3,1),sem);
        
        set(gca,'XTick',1:length(mB))
        set(gca,'XTickLabel',round(nanmedian(mB,1)*1)/1);
        xlabel('bearing')
        
        ylabel('reversal frequency (rev/s)')
        title (name)
        
    end
    leg{CC}=[num2str(round(delay*10/3)/10) 's'];
    CC=CC+1;
end

legend(h,leg);
%%
saveas(gca, 'reversal_freq.fig')





