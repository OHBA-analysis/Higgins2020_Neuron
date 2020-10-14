function fields_to_workspace(structure)

% fields_to_workspace(struct)
%
% adds all fields from the structure struct to the workspace
% will overwrite variables that already have corresponding field!
%
% returns field contents for 1x1 structure; cell array or matrix for
%   multi-dimensional structures (with same dimensions as structure)
%
% LH 081109/261012

if ~isstruct(structure)
    error('Input argument to fields_to_workspace must be structure');
end

s = size(structure);
names = fieldnames(structure);
if numel(structure)==1
    for n = 1:length(names)
        assignin('caller',names{n},getfield(structure,names{n}));
    end
else    %new section added 261012 for multidimensional structures:
   for n = 1:length(names)
       tmp = cell(s);
       for i = 1:numel(structure)
           [c{1:length(s)}] = ind2sub(s,i);
           tmp{i} = getfield(structure,c,names{n});
       end
       try
           tmp(cellfun(@isempty,tmp)) = {0};
           tmp = cell2mat(tmp);
       end
       assignin('caller',names{n},tmp);
   end
end