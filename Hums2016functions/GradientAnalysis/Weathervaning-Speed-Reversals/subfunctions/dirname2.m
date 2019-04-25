%outputs name of current directory +1 upper level as string
%
%(c) Julia Riedl, julia.riedl@imp.ac.at
%Created 2016

function name=dirname2(cd)

        nd=(cd);
        d= strfind(nd, '\');
        name=nd(d(end-1)+1:d(end)-1);
end