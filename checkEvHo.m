function checkEvHo(string)
%Litte function to visually inspect whether the 010219 tracking methods got
%accurate lawn boundaries.
if ~exist('string','var')
    string = '*010219*BACKGROUND.mat';
end
[~, pathname, ~] = uigetfile({'*'}); %this is nice because you can see which ones are in progress
cd(pathname);
D = rdir([pathname '\**\' string]);

for i = 1:length(D)
    bg_struct = load(D(i).name);
    bg_struct = bg_struct.bg_struct;
    figure();
    ydim = ceil(length(bg_struct)/3);
    [ha,~] = tight_subplot(ydim,3,0.01,0.01,0.01);
    set(gcf,'Position',[455 191 1200 900]);set(gcf,'Color','w');
    for j = 1:length(bg_struct)
        axes(ha(j)); hold on;
        imshow(bg_struct(j).orig_background);
        plot(bg_struct(j).ev_ho(:,1),bg_struct(j).ev_ho(:,2),'LineWidth',2);
        text(50,50,bg_struct(j).videoname,'FontSize',10,'Interpreter','none');
    end
end

end

