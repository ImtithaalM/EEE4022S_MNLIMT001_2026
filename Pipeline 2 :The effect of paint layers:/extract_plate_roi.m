function [roi, info] = extract_plate_roi(data, showPlot, label, quiet)
% EXTRACT_PLATE_ROI  Cuts the plate out of a thermal frame (outline A, B or C, whichever fits) and masks specks and cold pixels; large hot regions are kept.
% Outline A = original position, B = PAINT502 (IV_02306), C = PAINT1002 (IV_02216). Tried in order A, B, C: first with <= 5% background in an edge strip wins.
% Output: roi (NaN = masked) and info (rows, cols, position, sigma, hotAreaPct, pctCold, median, mean, std ...). showPlot = true draws the check figure.

%1.Outlines (rows = top, right, bottom, left; cols = x1 y1 x2 y2)
O.name = {'A','B','C'};
O.lines{1} = [56 59 579 59;  579 59 579 434;  56 434 579 434;  56 59 56 434];   % A (MATLAB 1-based)
O.lines{2} = [75 60.9 551 77.1;  540.6 424 546.2 74;  67 424 539 423.9;  69 419.9 74.1 58.9];   % B
O.lines{3} = [78 59 570 95;  560.8 439 565 93;  76 431.7 562 438;  72 434 81.9 58];            % C

%2.Settings
INSET_T = 28;  INSET_B = 15;  INSET_L = 35;  INSET_R = 20;   % px kept clear inside the outline

K_SIGMA = 5;  DILATE_PIX = 2;  SPECK_WIN = 31;  SPECK_FRAC = 0.10;   % outlier rule

stripPx = 6;  coldDelta = 4.5;  coldK = 2;  stripLimit = 0.2;  pickLimit = 0.05;  medDrop = 3;  medDropK = 1.5;   % outline pick: same strip test as check_roi_fit

%3.Pick outline
[nr, nc] = size(data);
nCand = numel(O.name);
score = inf(1, nCand);  ok = false(1, nCand);  R = cell(1, nCand);  C = cell(1, nCand);
for k = 1:nCand
    L = O.lines{k};
    r0 = ceil(max(L(1,[2 4]))) + INSET_T;   r1 = floor(min(L(3,[2 4]))) - INSET_B;
    c0 = ceil(max(L(4,[1 3]))) + INSET_L;   c1 = floor(min(L(2,[1 3]))) - INSET_R;
    if r0 < 1 || c0 < 1 || r1 > nr || c1 > nc || r1 - r0 < 2*stripPx || c1 - c0 < 2*stripPx, continue; end
    R{k} = r0:r1;  C{k} = c0:c1;
    sub = data(R{k}, C{k});
    med = median(sub(:), 'omitnan');
    sg  = 1.4826 * median(abs(sub(:) - med), 'omitnan');
    cutoff = med - max(coldDelta, coldK * sg);
    s = {sub(1:stripPx, :), sub(end-stripPx+1:end, :), sub(:, 1:stripPx), sub(:, end-stripPx+1:end)};
    frac = cellfun(@(x) mean(x(:) < cutoff), s);
    drop = cellfun(@(x) med - median(x(:), 'omitnan'), s);
    score(k) = max(frac);
    ok(k) = max(frac) <= stripLimit && max(drop) <= max(medDrop, medDropK * sg);
end
% first clean outline; else the passing one with least background; else the lowest overall (check_roi_fit flags it)
pick = find(score <= pickLimit & ok, 1, 'first');
if isempty(pick)
    sc = score;  sc(~ok) = inf;
    [m, pick] = min(sc);
    if isinf(m), [~, pick] = min(score); end
end
if isinf(score(pick)), error('extract_plate_roi:noOutline', 'No outline fits inside the frame.'); end
rows = R{pick};  cols = C{pick};
roi  = data(rows, cols);

%4.Mask (median + MAD)
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

v = roi(~isnan(roi));
info.rows        = rows;
info.cols        = cols;
info.position    = O.name{pick};
info.scratchMask = bad;
info.sigma       = sigma;
info.pctHot      = 100 * nnz(hot)  / numel(hot);
info.hotAreaPct  = 100 * nnz(hot & ~speck) / numel(hot);
info.pctCold     = 100 * nnz(cold) / numel(cold);
info.pctFlagged  = pctFlagged;
info.pctMasked   = 100 * nnz(bad) / numel(bad);
info.mean        = mean(v);
info.median      = median(v);
info.std         = std(v);
info.max         = max(v);
info.min         = min(v);

if nargin < 3 || isempty(label), label = 'frame'; end
if nargin < 4 || isempty(quiet), quiet = false; end
if ~quiet
fprintf('\n=== %s ===\n', label);
fprintf('ROI    : outline %s, rows %d-%d, cols %d-%d  (%d x %d)\n', info.position, ...
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

%5.Check figure
if nargin > 1 && showPlot
    L = O.lines{pick};
    figure('Name', 'extract_plate_roi check');
    subplot(1,2,1);
    imagesc(data, [20 45]); axis image; colormap jet; colorbar; hold on;
    plot([L(1,1) L(1,3) L(3,3) L(3,1) L(1,1)], [L(1,2) L(1,4) L(3,4) L(3,2) L(1,2)], 'w-', 'LineWidth', 1.5);
    rectangle('Position', [cols(1) rows(1) numel(cols) numel(rows)], ...
        'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');
    title(sprintf('Outline %s (white), dashed = ROI used', info.position));
    subplot(1,2,2);
    imagesc(bad); axis image; colormap(gca, gray);
    title(sprintf('Masked specks/cold (%.2f%% of ROI)', info.pctMasked));
end
end