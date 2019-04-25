function changefilenames_recursively( topdir, toreplace_string, replacewith_string )
%little script to change the filenames in a directory based on comparing a
%toreplace_string and a replacewith_string. 05/09/2018

D = rdir([topdir '\**\*' toreplace_string '*']);

for i = 1:length(D)
    filename = D(i).name;
    startInd = strfind(filename,toreplace_string);
    endInd = startInd+length(toreplace_string)-1;
    for j = 1:size(startInd,2)
        filename(startInd(j):endInd(j)) = replacewith_string;
    end
    movefile(D(i).name,filename);
end

end

