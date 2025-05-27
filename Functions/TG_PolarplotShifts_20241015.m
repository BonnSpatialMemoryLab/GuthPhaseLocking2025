function TG_PolarplotShifts_20241015(alpha, beta, colorAlpha, colorBeta, colorLine, param)
% TG_POLARPLOTSHIFTS_20241015 creates a polarplot for angular shifts

circDiff            = angdiff(alpha, beta);

% sort by shift magnitude
[~, sortIdx]        = sort(abs(circDiff));
sortIdx             = flip(sortIdx);
circDiff            = circDiff(sortIdx);
alpha               = alpha(sortIdx);
beta                = beta(sortIdx);

% loop through all shifting cells
for iShift = 1:size(circDiff, 1)

    % plot angular difference
    shiftLine       = linspace(alpha(iShift), alpha(iShift) + circDiff(iShift), 200);
    shiftRho        = repmat(iShift, size(shiftLine, 2), 1);
    polarplot(shiftLine, shiftRho, 'color', colorLine, 'LineStyle', '-', 'LineWidth', 1);
    hold on;
    set(gca, 'ThetaTickLabel', param.polarHistLabels);
    polarscatter(alpha(iShift), iShift, '.', 'MarkerEdgeColor', colorAlpha, 'SizeData', 200);
    polarscatter(beta(iShift), iShift, '.', 'MarkerEdgeColor', colorBeta, 'SizeData', 200);
end
rticks([]);
rlim([0, size(circDiff, 1)]);

end