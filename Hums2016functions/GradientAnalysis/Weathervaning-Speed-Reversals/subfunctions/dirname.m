%outputs name of current directory as string
%
%(c) Julia Riedl, julia.riedl@imp.ac.at
%Created 2016

function name=dirname(cd)

        nd=(cd);
        d= strfind(nd, '\');
        name=nd(d(end)+1:end);
end