function [] = f_save_mat_fig(h,filename)
saveas(h,filename)
map = load(filename,'-mat');
names = fieldnames(map);
for j = 1:numel(names)
    map.(names{j}).properties.Visible = 'on';
end
save(filename,'-struct','map');
end