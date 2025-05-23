function unwrappedCell = tg_unwrapNestedCell(inNestedCell)
%
% tg_unwrapNestedCell unwraps a nested cell array.

unwrappedCell   = cell(sum(cellfun(@(x) numel(x), inNestedCell)), 1);
idx             = 1;
for i = 1:size(inNestedCell, 1)
    for j = 1:numel(inNestedCell{i})
        unwrappedCell{idx, 1}   = inNestedCell{i, 1}{j};
        idx                     = idx + 1;
    end
end