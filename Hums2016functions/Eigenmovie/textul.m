function texthandle=textul(txt,fontsize,color)
%textul(txt,fontsize)
%
% write text in the upper left corner of a plot
%
% convenience function
% (c) Saul Kato, saul@kato.com 10-07-10


if nargin<3
    color='k';
end

if nargin<2
    fontsize=10;
end

ax=axis;
texthandle=text(ax(1),ax(4),txt,...
            'VerticalAlignment','top','HorizontalAlignment','left',...
            'FontSize',fontsize,'Color',color);