function ppc = TG_PPC_20241128(phases)

% TG_PPC_20241128 computes pairwise phase consistency (PPC)
%
% Input:
%   phases    vector of phases
%
% Output:
%   ppc       pairwise phase consistency
%
% References:
%   Vinck et al., NeuroImage, 2010
%   https://gist.github.com/mrkrause/cde35bc5a5cdbb806773e914dce186b3
%
% Tim Guth, 2024

% number of phases
N               = length(phases);

% if fewer than two phases
if N < 2
    ppc             = NaN;
    warning('Not enough phases to compute pairwise phase consistency');
    return;
end

% scale factor
scaleFactor     = 2 / (N * (N - 1));

% sine and cosine of phases
sinPhases       = sin(phases);
cosPhases       = cos(phases);

% compute pairwise phase consistency
ppc             = 0;
for iPhase = 1:(N - 1)

    % pairwise cosine-sine product summation
    pairwiseSum     = sum(...
        cosPhases(iPhase) .* cosPhases((iPhase + 1):N) + ...
        sinPhases(iPhase) .* sinPhases((iPhase + 1):N) ...
        );

    % increment PPC
    ppc             = ppc + scaleFactor * pairwiseSum;
end

end