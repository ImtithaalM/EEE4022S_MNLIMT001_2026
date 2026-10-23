function [roi, info] = extract_depth_roi(data, showPlot, label, quiet, plate)
% EXTRACT_DEPTH_ROI  Cuts the plate out of a thermal frame and measures the depth rectangles of depth_geometry(plate) ('2mm' default, '4mm').
% roi = plate pixels (NaN outside the plate and at hot specks). info = rows, cols, plateMask, mask, dist, refMedian, refSigma, rect(k) (median, dT, pctValid ...), rectFit.
% showPlot = true draws the check figure, quiet = true suppresses the printout. Cold pixels are kept: bare steel reads colder and that is the signal.

if nargin < 2 || isempty(showPlot), showPlot = false; end
if nargin < 3 || isempty(label),    label    = 'frame'; end
if nargin < 4 || isempty(quiet),    quiet    = false; end
if nargin < 5 || isempty(plate),    plate    = '2mm'; end

g  = depth_geometry(plate);
[nr, nc] = size(data);
nR = size(g.rectBox, 1);

%1.Distance inside each outline line (+ = inside the plate)
[X, Y] = meshgrid(1:nc, 1:nr);
edges = [g.top; g.right; g.bottom; g.left];
cx = mean((edges(:,1) + edges(:,3)) / 2);        % a point well inside the plate
cy = mean((edges(:,2) + edges(:,4)) / 2);
dist = zeros(nr, nc, 4);
for k = 1:4
    dx = edges(k,3) - edges(k,1);
    dy = edges(k,4) - edges(k,2);
    s  = sign(dx * (cy - edges(k,2)) - dy * (cx - edges(k,1)));   % which side the plate is on
    dist(:,:,k) = s * (dx * (Y - edges(k,2)) - dy * (X - edges(k,1))) / hypot(dx, dy);
end

%2.Plate ROI (inside all four lines by the inset)
plateMask = true(nr, nc);
for k = 1:4
    plateMask = plateMask & (dist(:,:,k) >= g.inset(k));
end
rr = find(any(plateMask, 2));
cc = find(any(plateMask, 1));
rows = rr(1):rr(end);
cols = cc(1):cc(end);

%3.Hot-speck mask
sub    = data(rows, cols);
inMask = plateMask(rows, cols);
vals   = sub(inMask);
med    = median(vals, 'omitnan');
sigma  = 1.4826 * median(abs(vals - med), 'omitnan');
hot    = inMask & (sub > med + g.kSigma * sigma);
box    = ones(g.speckWin, 1) / g.speckWin;
hotDensity = conv2(box, box', double(hot), 'same');
speck  = hot & (hotDensity < g.speckFrac);
badSub = conv2(double(speck), ones(2*g.dilatePix + 1), 'same') > 0;
bad    = false(nr, nc);
bad(rows, cols) = badSub;
validMask = plateMask & ~bad;

roi = sub;
roi(~inMask | badSub) = NaN;

%4.Rectangle masks and reference region
rectMask = false(nr, nc, nR);
ext = false(nr, nc);
for k = 1:nR
    b  = g.rectBox(k, :);   ri = g.rectInset;
    rectMask(b(3)+ri:b(4)-ri, b(1)+ri:b(2)-ri, k) = true;
    ext(max(1, b(3)-g.refMargin):min(nr, b(4)+g.refMargin), ...
        max(1, b(1)-g.refMargin):min(nc, b(2)+g.refMargin)) = true;
end
for k = 1:size(g.refExclude, 1)        % bare-steel patches left out of the reference
    b = g.refExclude(k, :);
    ext(max(1, b(3)):min(nr, b(4)), max(1, b(1)):min(nc, b(2))) = true;
end
refMask   = validMask & ~ext;
refVals   = data(refMask);
refMedian = median(refVals, 'omitnan');
refSigma  = 1.4826 * median(abs(refVals - refMedian), 'omitnan');

%5.Rectangle statistics
rect    = struct('name', {}, 'nTotal', {}, 'nValid', {}, 'pctValid', {}, 'mean', {}, 'median', {}, ...
                 'std', {}, 'min', {}, 'max', {}, 'dT', {}, 'dTmean', {});
rectFit = zeros(1, nR);
for k = 1:nR
    m = rectMask(:,:,k);
    v = data(m & validMask);
    nValid = numel(v);
    if nValid == 0, v = NaN; end
    rect(k).name     = g.rectNames{k};
    rect(k).nTotal   = nnz(m);
    rect(k).nValid   = nValid;
    rect(k).pctValid = 100 * nValid / nnz(m);
    rect(k).mean     = mean(v);
    rect(k).median   = median(v);
    rect(k).std      = std(v);
    rect(k).min      = min(v);
    rect(k).max      = max(v);
    rect(k).dT       = rect(k).median - refMedian;
    rect(k).dTmean   = rect(k).mean   - refMedian;

    % clearance to the plate ROI edge
    b = g.rectBox(k, :);
    idx = sub2ind([nr nc], [b(3) b(3) b(4) b(4)], [b(1) b(2) b(1) b(2)]);   % the 4 corners
    clear_ = inf;
    for kk = 1:4
        dk = dist(:,:,kk);
        clear_ = min(clear_, min(dk(idx)) - g.inset(kk));
    end
    rectFit(k) = clear_;
end

%6.Whole-plate statistics
v = data(validMask);
info.rows        = rows;
info.cols        = cols;
info.plateMask   = plateMask;
info.mask        = validMask;
info.scratchMask = bad;
info.dist        = dist;
info.inset       = g.inset;
info.rectMask    = rectMask;
info.refMedian   = refMedian;
info.refSigma    = refSigma;
info.mean        = mean(v);
info.median      = median(v);
info.std         = std(v);
info.max         = max(v);
info.min         = min(v);
info.pctMasked   = 100 * nnz(bad & plateMask) / nnz(plateMask);
info.rectNames   = g.rectNames;
info.rect        = rect;
info.rectFit     = rectFit;

if ~quiet
    fprintf('\n=== %s ===\n', label);
    fprintf('Plate ROI : rows %d-%d, cols %d-%d (bounding box), %d px, %.2f%% hot specks masked\n', ...
        rows(1), rows(end), cols(1), cols(end), nnz(plateMask), info.pctMasked);
    fprintf('Reference : median %.2f | sigma %.2f degC  (plate outside the rectangles)\n', refMedian, refSigma);
    fprintf('%-5s %8s %8s %8s %8s %8s %7s %9s\n', 'Rect', 'Median', 'Mean', 'Std', 'dT', 'dTmean', 'Valid%', 'Fit(px)');
    for k = 1:nR
        fprintf('%-5s %8.2f %8.2f %8.2f %8.2f %8.2f %7.1f %9.1f\n', rect(k).name, rect(k).median, ...
            rect(k).mean, rect(k).std, rect(k).dT, rect(k).dTmean, rect(k).pctValid, rectFit(k));
    end
end

%7.Check figure
if showPlot
    figure('Name', ['extract_depth_roi check: ' label]);
    subplot(1,2,1);
    imagesc(data, [refMedian - 30, refMedian + 8]); axis image; colormap(gca, jet); colorbar; hold on;
    for k = 1:4
        plot([edges(k,1) edges(k,3)], [edges(k,2) edges(k,4)], 'w-', 'LineWidth', 1.5);
    end
    contour(double(plateMask), [0.5 0.5], 'k--', 'LineWidth', 1.5);
    for k = 1:nR
        b = g.rectBox(k, :);
        rectangle('Position', [b(1) b(3) b(2)-b(1) b(4)-b(3)], 'EdgeColor', 'm', 'LineWidth', 1.5);
        text(b(1) + 4, b(3) + 12, g.rectNames{k}, 'Color', 'm', 'FontWeight', 'bold');
    end
    title('White = plate outline, dashed = plate ROI, magenta = rectangles');
    subplot(1,2,2);
    lab = zeros(nr, nc);
    lab(plateMask) = 1;
    for k = 1:nR, lab(rectMask(:,:,k)) = 1 + k; end
    lab(bad & plateMask) = 6;
    imagesc(lab, [0 6]); axis image;
    colormap(gca, [0.15 0.15 0.15; 0.65 0.65 0.65; 0.2 0.6 1; 1 0.5 0.1; 0.2 0.8 0.3; 0.9 0.2 0.8; 1 1 0]);
    title(sprintf('Grey = reference, colours = rectangles, yellow = masked specks (%.2f%%)', info.pctMasked));
end
end