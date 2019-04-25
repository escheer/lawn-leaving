function copyfiles_recursively_mixedRegen( topdir, outpath, string )
%little script to copy all files matching a regular expression to an
%outpath

D = rdir([topdir '\**\' string]);
filenames = {D(:).name}';
regenfiles = contains(filenames,'REGEN');
filenames(find(regenfiles)-1) = [];

for i = 1:length(filenames)
    copyfile(filenames{i},outpath);
end

end

