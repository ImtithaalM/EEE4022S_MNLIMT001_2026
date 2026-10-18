function [roi, info] = extract_plate_roi(data, showPlot, label, quiet)
% EXTRACT_PLATE_ROI  Cuts the plate out of a 480x640 thermal frame and masks glints, pits and cold pixels (masked = NaN). Large hot regions are kept.
% [roi, info] = extract_plate_roi(data, showPlot, label, quiet): showPlot = draw check figure, label = printout name, quiet = no printout.
% The outline is hard-coded from IRSoft line profiles (camera and plate fixed); redo it only if either moves.

% ---- plate outline ----
OUT_TOP    = 59;
OUT_BOTTOM = 434;
OUT_LEFT   = 56;
OUT_RIGHT  = 579;

% ---- inset (corrosion reaching the plate edge is partly cropped) ----
INSET_T = 28;
INSET_B = 15;
INSET_L = 35;
INSET_R = 20;

% ---- outlier rule ----
K_SIGMA     = 5;      % flag pixels more than K * robust sigma from the median
DILATE_PIX  = 2;      % grow each masked pixel by this many pixels
SPECK_WIN   = 31;     % window (px) used to tell specks from patches
SPECK_FRAC  = 0.10;   % a hot pixel is a speck if fewer than this fraction of its window is hot

rows = (OUT_TOP + INSET_T) : (OUT_BOTTOM - INSET_B);
cols = (OUT_LEFT + INSET_L) : (OUT_RIGHT - INSET_R);
roi  = data(rows, cols);

%1.Mask
med   = median(roi(:), 'omitnan');
sigma = 1.4826 * median(abs(roi(:) - med), 'omitnan');
hot   = roi > med + K_SIGMA * sigma;
cold  = roi < med - K_SIGMA * sigma;
box   = ones(SPECK_WIN, 1) / SPECK_WIN;
hotDensity = conv2(box, box', double(hot), 'same');
speck = hot & (hotDensity < SPECK_FRAC);
bad   = speck | cold;
pctFlagged = 100 * nnz(bad) / numel(bad);
kern  = ones(2*DILATE_PIX + 1);
bad   = conv2(double(bad), kern, 'same') > 0;
roi(bad) = NaN;

%2.Stats
v = roi(~isnan(roi));
info.rows        = rows;
info.cols        = cols;
info.scratchMask = bad;
info.sigma       = sigma;                                  % robust spread of the raw ROI (not inflated by corrosion)
info.pctHot      = 100 * nnz(hot)  / numel(hot);
info.hotAreaPct  = 100 * nnz(hot & ~speck) / numel(hot);   % % of ROI in large hot regions
info.pctCold     = 100 * nnz(cold) / numel(cold);
info.pctFlagged  = pctFlagged;
info.pctMasked   = 100 * nnz(bad) / numel(bad);
info.mean        = mean(v);
info.median      = median(v);
info.std         = std(v);
info.max         = max(v);
info.min         = min(v);

%3.Printout
if nargin < 3 || isempty(label), label = 'frame'; end
if nargin < 4 || isempty(quiet), quiet = false; end
if ~quiet
fprintf('\n=== %s ===\n', label);
fprintf('ROI    : rows %d-%d, cols %d-%d  (%d x %d)\n', ...
    rows(1), rows(end), cols(1), cols(end), numel(rows), numel(cols));
fprintf('Masked : %.2f%% specks/cold flagged -> %.2f%% after dilation; %.2f%% large hot region kept\n', ...
    info.pctFlagged, info.pctMasked, info.hotAreaPct);
fprintf('Plate  : median %.2f | mean %.2f | std %.2f | range %.2f to %.2f degC\n', ...
    info.median, info.mean, info.std, info.min, info.max);
end

if info.pctFlagged > 3 && ~quiet
    warning('extract_plate_roi:manyFlagged', ...
        'Over 3%% of pixels flagged as specks/cold - check K_SIGMA or look at the frame.');
end

%4.Check figure
if nargin > 1 && showPlot
    figure('Name', 'extract_plate_roi check');
    subplot(1,2,1);
    imagesc(data, [20 45]); axis image; colormap jet; colorbar; hold on;
    rectangle('Position', [OUT_LEFT OUT_TOP OUT_RIGHT-OUT_LEFT OUT_BOTTOM-OUT_TOP], ...
        'EdgeColor', 'w', 'LineWidth', 1.5);
    rectangle('Position', [cols(1) rows(1) numel(cols) numel(rows)], ...
        'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');
    title('White = plate outline, dashed = ROI used');
    subplot(1,2,2);
    imagesc(bad); axis image; colormap(gca, gray);
    title(sprintf('Masked specks/cold (%.2f%% of ROI)', info.pctMasked));
end
end
