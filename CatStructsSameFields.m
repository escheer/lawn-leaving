% function Z = CatStructsSameFields (S,T)
% 
% Z = cell2struct(cellfun(@vertcat,struct2cell(S),struct2cell(T),'uni',0),fieldnames(S),1);
% 
% end

function S = CatStructsSameFields(S, T, dim)
fields = fieldnames(S);
for k = 1:numel(fields)
%   disp(k);
%   if k == 21
%       disp('debug');
%   end
  aField     = fields{k};
  S.(aField) = cat(dim, S.(aField), T.(aField));
end