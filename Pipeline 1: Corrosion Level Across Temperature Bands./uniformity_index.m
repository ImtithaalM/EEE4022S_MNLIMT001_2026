function idx = uniformity_index(stdRatio, stdRatioAnchor)
% UNIFORMITY_INDEX  Linear rescale of the std ratio: 0 = like bare steel (ratio ~1), 1 = like 100% corroded (ratio = anchor). Not clipped.
% stdRatio = test spread / bare spread (compute_indicators_multi); stdRatioAnchor = same ratio for fully corroded frames of the same session.
% A 2-point score calibrated only at 0% and 100%, not a corrosion percentage.
idx = (1 - stdRatio) ./ (1 - stdRatioAnchor);
end
