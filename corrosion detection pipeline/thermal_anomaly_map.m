% THERMAL_ANOMALY_MAP  Thermal-only anomaly map and area % for bare (control), light, FULL and HALF on the same plate ROI, per band group.
% Removes the smooth plate pattern (two-pass Gaussian background), flags pixels beyond kThresh robust sigmas of the same-band bare frames, and also builds a steel-referenced hot map (use that one for HALF).
% Needs frame_bands_<plate>.csv from check_roi_fit. Output: stats\anomaly_results_withTape.csv, anomaly_masks_withTape.mat and one figure per band group.

% ---- settings ----
%root = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\matlab pipeline';
root        = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\matlab pipeline';
statsFolder = fullfile(root, 'stats');
testPlates  = ["light" "FULL" "HALF"];
% tag, band labels (must match the CSVs), bare control folder / CSV / frames kept
bandGroups = struct( ...
    'tag',      {'22-24', '24-26', '28-32', '38-42', '48-52'}, ...
    'bands',    {"22-24", "24-26", ["28-30" "30-32"], ["38-40" "40-42"], "48-52"}, ...
    'bareDir',  {fullfile(root, 'BARE NEW'), fullfile(root, 'BARE NEW'), fullfile(root, 'BARE NEW'), fullfile(root, 'BARE NEW'), fullfile(root, 'BARE NEW')}, ...
    'bareCsv',  {'frame_bands_BARE NEW.csv', 'frame_bands_BARE NEW.csv', 'frame_bands_BARE NEW.csv', 'frame_bands_BARE NEW.csv', 'frame_bands_BARE NEW.csv'}, ...
    'bareKeep', {2172:2192, 2172:2192, 2142:2153, 2115:2123, 2097:2102});
useGroups   = {'22-24', '24-26', '28-32', '38-42', '48-52'};   % groups to run
mergeBands  = ["48-50" "50-52"];   % CSV bands merged into one ...
mergedLabel = "48-52";             % ... with this label
excludeFrames = struct('HALF', 1573);   % test-plate frames to drop
minBareFrames = 2;            % a band needs at least this many bare frames
sigmaBG     = 30;             % Gaussian scale (px) of the removed background
bgTwoPass   = true;           % leave the warm region out of the background
bgGrow      = 5;              % px the warm region is grown by before it is left out
kThresh     = 5;              % anomaly = residual beyond kThresh * sigma(bare)
minBlob     = 5;              % ignore blobs smaller than this (px)
refSmooth   = 2;              % px smoothing for the steel-referenced map
refLowPct   = 20;  refHighPct = 95;   % steel level percentile, and percentile that sets DELTA on bare
closeRad    = 4;              % px: warm pixels closer than this are joined
bigBlob     = 500;            % px: joined region at least this big = large region
excludeTape   = false;        % set the masking-tape area to NaN
excludeCorner = false;        % also exclude the hot bottom-left corner
cornerPoly  = [5 272; 110 278; 110 335; 5 335];   % [column row] corners of that corner
% frames drawn per band group: group, plate, frame (falls back to the plate's first frame in the group)
exampleFrames = {
    '22-24', 'light', 1898;  '22-24', 'HALF', 1719;  '22-24', 'bare (control)', 2185;
    '24-26', 'light', 1883;  '24-26', 'FULL', 1770;  '24-26', 'HALF', 1715;  '24-26', 'bare (control)', 2174;
    '28-32', 'light', 1853;  '28-32', 'FULL', 1756;  '28-32', 'HALF', 1680;  '28-32', 'bare (control)', 2148;
    '38-42', 'light', 1810;  '38-42', 'FULL', 1749;  '38-42', 'HALF', 1620;  '38-42', 'bare (control)', 2120;
    '48-52', 'light', 1786;  '48-52', 'FULL', 1737;  '48-52', 'HALF', 1580;  '48-52', 'bare (control)', 2099};
tag = '';
if ~excludeTape,  tag = [tag '_withTape'];  end
if excludeCorner, tag = [tag '_noCorner'];  end

%1.Exclusion mask
extraMask = [];
if excludeTape || excludeCorner
    S = load(fullfile(statsFolder, 'partial_patch_masks.mat'), 'tapeMask', 'refSize');
    extraMask = false(S.refSize);
    if excludeTape,   extraMask = extraMask | S.tapeMask; end
    if excludeCorner, extraMask = extraMask | poly2mask(cornerPoly(:,1), cornerPoly(:,2), S.refSize(1), S.refSize(2)); end
end

%2.Frame lists
jobs = cell(0, 5);      % plate label, folder, frame, band, band-group tag
targetBands = strings(1, 0);
for g = 1:numel(bandGroups)
    if ~ismember(bandGroups(g).tag, useGroups), continue; end
    gb = bandGroups(g).bands;  targetBands = [targetBands, gb]; %#ok<AGROW>
    for p = 1:numel(testPlates)
        Tp = readtable(fullfile(statsFolder, ['frame_bands_' char(testPlates(p)) '.csv']), 'TextType', 'string');
        Tp.band(ismember(strtrim(Tp.band), mergeBands)) = mergedLabel;
        Tp = splitBand(Tp);
        ok = startsWith(string(Tp.status), "OK") & ismember(strtrim(Tp.band), gb);
        if isfield(excludeFrames, char(testPlates(p)))
            drop = ok & ismember(Tp.frame, excludeFrames.(char(testPlates(p))));
            if any(drop), fprintf('%s group %s: excluding frame(s) %s\n', testPlates(p), bandGroups(g).tag, mat2str(Tp.frame(drop)')); end
            ok = ok & ~drop;
        end
        fr = Tp.frame(ok);  bd = strtrim(Tp.band(ok));
        for i = 1:numel(fr)
            jobs(end+1, :) = {char(testPlates(p)), fullfile(root, char(testPlates(p))), fr(i), char(bd(i)), bandGroups(g).tag}; %#ok<SAGROW>
        end
    end
    Tb = readtable(fullfile(statsFolder, bandGroups(g).bareCsv), 'TextType', 'string');
    Tb.band(ismember(strtrim(Tb.band), mergeBands)) = mergedLabel;
    Tb = splitBand(Tb);
    okB = startsWith(string(Tb.status), "OK") & ismember(strtrim(Tb.band), gb);
    if ~isempty(bandGroups(g).bareKeep), okB = okB & ismember(Tb.frame, bandGroups(g).bareKeep); end
    frB = Tb.frame(okB);  bdB = strtrim(Tb.band(okB));
    for i = 1:numel(frB)
        jobs(end+1, :) = {'bare (control)', bandGroups(g).bareDir, frB(i), char(bdB(i)), bandGroups(g).tag}; %#ok<SAGROW>
    end
end
keepJob = true(size(jobs, 1), 1);
for j = 1:size(jobs, 1)
    nb = nnz(strcmp(jobs(:,1), 'bare (control)') & strcmp(jobs(:,4), jobs{j,4}));
    if nb < minBareFrames, keepJob(j) = false; end
end
for bnd = targetBands
    nb = nnz(strcmp(jobs(:,1), 'bare (control)') & strcmp(jobs(:,4), char(bnd)));
    if nb < minBareFrames, fprintf('Band %s skipped: %d bare frame(s), need %d.\n', bnd, nb, minBareFrames); end
end
jobs = jobs(keepJob, :);
n = size(jobs, 1);
if n == 0, error('No frames to process: check bandGroups against the band labels in the CSVs.'); end
fprintf('%d frames to process (%d bare control).\n', n, nnz(strcmp(jobs(:,1), 'bare (control)')));

%3.Residual map of every frame
res = cell(n, 1);  plateMed = nan(n, 1);
tsm = cell(n, 1);  pLow = nan(n, 1);  pHigh = nan(n, 1);
for j = 1:n
    roi = loadROI(jobs{j,2}, jobs{j,3});
    if ~isempty(extraMask)
        if ~isequal(size(roi), size(extraMask)), error('Frame %d ROI size differs from the mask - redraw the masks.', jobs{j,3}); end
        roi(extraMask) = NaN;
    end
    plateMed(j) = median(roi(:), 'omitnan');
    res{j} = highpass(roi, sigmaBG, bgTwoPass, bgGrow);
    vv = ~isnan(roi);  zz = roi;  zz(~vv) = median(roi(vv));
    ts = imgaussfilt(zz, refSmooth);  ts(~vv) = NaN;
    tsm{j} = single(ts);  pLow(j) = prctile(ts(vv), refLowPct);  pHigh(j) = prctile(ts(vv), refHighPct);
end

%4.Sigma per frame (same-band bare residuals, leave-one-out for bare)
isBare = strcmp(jobs(:,1), 'bare (control)');
sig = nan(n, 1);
for j = 1:n
    ref = find(isBare & strcmp(jobs(:,4), jobs{j,4}) & (1:n)' ~= j);
    if isempty(ref), continue; end
    pooled = cell2mat(cellfun(@(x) x(~isnan(x)), res(ref), 'UniformOutput', false));
    sig(j) = 1.4826 * median(abs(pooled - median(pooled)));
end

%5.Steel-referenced DELTA (one per band group, never pooled)
groupTags = unique(jobs(:,5), 'stable');
deltaRefG = zeros(numel(groupTags), 1);
for g = 1:numel(groupTags)
    mg = isBare & strcmp(jobs(:,5), groupTags{g});
    deltaRefG(g) = median(pHigh(mg) - pLow(mg));
    fprintf('Steel-referenced map, band group %s: DELTA = %.2f degC above the steel level (%dth percentile), set from %d bare frames.\n', groupTags{g}, deltaRefG(g), refLowPct, nnz(mg));
end

%6.Anomaly masks and area %
rows = cell(n, 1);  hotM = cell(n, 1);  coldM = cell(n, 1);  hotRefM = cell(n, 1);
for j = 1:n
    valid = ~isnan(res{j});
    hotRefM{j} = bwareaopen(double(tsm{j}) > pLow(j) + deltaRefG(strcmp(groupTags, jobs{j,5})), minBlob);
    hotM{j}  = bwareaopen(res{j} >  kThresh * sig(j), minBlob);
    coldM{j} = bwareaopen(res{j} < -kThresh * sig(j), minBlob);
    nv = nnz(valid);
    pk = NaN;  if any(hotM{j}(:)), pk = max(res{j}(hotM{j})); end
    cc = bwconncomp(imclose(hotM{j}, strel('disk', closeRad)));
    bsz = cellfun(@numel, cc.PixelIdxList);
    largestPct = 0;  if ~isempty(bsz), largestPct = 100 * max(bsz) / nv; end
    bigPct = 100 * sum(bsz(bsz >= bigBlob)) / nv;
    nBig   = nnz(bsz >= bigBlob);
    % edge-pixel-to-area ratio (perimeter pixels / region pixels)
    edgeAll = NaN;  if any(hotM{j}(:)), edgeAll = nnz(bwperim(hotM{j})) / nnz(hotM{j}); end
    bigMask = false(size(hotM{j}));  idxBig = find(bsz >= bigBlob);
    for q = idxBig(:)', bigMask(cc.PixelIdxList{q}) = true; end
    edgeBig = NaN;  if any(bigMask(:)), edgeBig = nnz(bwperim(bigMask)) / nnz(bigMask); end
    edgeRef = NaN;  if any(hotRefM{j}(:)), edgeRef = nnz(bwperim(hotRefM{j})) / nnz(hotRefM{j}); end
    rows{j} = {jobs{j,1}, jobs{j,4}, jobs{j,3}, plateMed(j), sig(j), ...
               100 * nnz(hotM{j} | coldM{j}) / nv, 100 * nnz(hotM{j}) / nv, 100 * nnz(coldM{j}) / nv, pk, largestPct, bigPct, nBig, 100 * nnz(hotRefM{j}) / nv, edgeAll, edgeBig, edgeRef};
end
R = cell2table(vertcat(rows{:}), 'VariableNames', {'Plate','Band','Frame','PlateMedian','SigmaBare','AnomPct','HotPct','ColdPct','PeakC','LargestPct','BigBlobPct','NBig','HotRefPct','EdgeAreaAll','EdgeAreaBig','EdgeAreaRef'});

%7.Save
writetable(R, fullfile(statsFolder, ['anomaly_results' tag '.csv']));
resS = cellfun(@single, res, 'UniformOutput', false);
save(fullfile(statsFolder, ['anomaly_masks' tag '.mat']), 'jobs', 'resS', 'hotM', 'coldM', 'hotRefM', 'groupTags', 'deltaRefG', 'sig', 'R');

%8.Summary
fprintf('\nBackground scale %d px (two-pass %d), threshold %.1f sigma(bare), blobs >= %d px. Tape excluded: %d, corner excluded: %d.\n', ...
        sigmaBG, bgTwoPass, kThresh, minBlob, excludeTape, excludeCorner);
fprintf('%-15s %-8s %-3s %-22s %-8s %-8s %-9s %-7s %-11s %-9s %-10s %-9s\n', 'Plate', 'Band', 'N', '<strong>ANOMALY AREA %</strong>', 'Hot %', 'Cold %', 'SigmaBare', 'PeakC', 'LargestPct', 'BigBlob%', 'SteelRef%', 'Edge/Area');
for p = unique(R.Plate, 'stable')'
    for bnd = targetBands
        m = strcmp(R.Plate, p{1}) & strcmp(R.Band, bnd);
        if ~any(m), continue; end
        fprintf('%-15s %-8s %-3d <strong>%6.2f</strong> (sd %5.2f)     %-8.2f %-8.2f %-9.3f %-7.1f %-11.2f %-9.2f %-10.1f %-9.3f\n', p{1}, bnd, nnz(m), ...
            mean(R.AnomPct(m)), std(R.AnomPct(m)), mean(R.HotPct(m)), mean(R.ColdPct(m)), mean(R.SigmaBare(m)), mean(R.PeakC(m), 'omitnan'), mean(R.LargestPct(m)), mean(R.BigBlobPct(m)), mean(R.HotRefPct(m)), mean(R.EdgeAreaAll(m), 'omitnan'));
    end
end
fprintf('\n');
disp(R);

%9.Figures (one per band group)
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
for tg = groupTags'
    ex = exampleFrames(strcmp(exampleFrames(:,1), tg{1}), 2:3);
    f = figure('Name', ['Anomaly maps ' tg{1}], 'Position', [30 30 1500 1000]);
    for e = 1:size(ex, 1)
        inG = strcmp(jobs(:,1), ex{e,1}) & strcmp(jobs(:,5), tg{1});
        j = find(inG & cell2mat(jobs(:,3)) == ex{e,2}, 1);
        if isempty(j), j = find(inG, 1); end
        if isempty(j), continue; end
        ex{e,2} = jobs{j,3};
        v = ~isnan(res{j});
        subplot(4, 3, 3*e - 2);
        h = imagesc(res{j}, [-5 5] * sig(j)); set(h, 'AlphaData', v); set(gca, 'Color', [0.85 0.85 0.85]); axis image; colorbar;
        title(sprintf('%s IV\\_%05d: local residual (degC)', ex{e,1}, ex{e,2}));
        subplot(4, 3, 3*e - 1);
        h = imagesc(double(hotM{j}) - double(coldM{j}), [-1 1]); set(h, 'AlphaData', v); set(gca, 'Color', [0.85 0.85 0.85]); axis image;
        title(sprintf('Anomalies %.1f%% (yellow hot, dark blue cold)', R.AnomPct(j)));
        subplot(4, 3, 3*e);
        h = imagesc(double(hotRefM{j}), [0 1]); set(h, 'AlphaData', v); set(gca, 'Color', [0.85 0.85 0.85]); axis image;
        title(sprintf('Steel-referenced hot area %.1f%% (yellow = hotter than steel + DELTA)', R.HotRefPct(j)));
    end
    try, exportgraphics(f, fullfile(statsFolder, ['anomaly_maps_' tg{1} tag '_' stamp '.png']), 'Resolution', 200); catch, end
end

% ---------------------------------------------------------------------
function T = splitBand(T)
% Splits rows labelled '22-26' into '22-24' (plate median 22 to <24) and '24-26' (24 to <26); other labels are left alone.
b = strtrim(string(T.band));
T.band = b;
T.band(b == "22-26" & T.median >= 22 & T.median < 24) = "22-24";
T.band(b == "22-26" & T.median >= 24 & T.median < 26) = "24-26";
end

function r = highpass(roi, s, twoPass, grow)
% Local residual = pixel minus a NaN-aware Gaussian background of scale s px; twoPass leaves the warm region out of the background.
valid = ~isnan(roi);
[bg1, ~] = blurbg(roi, s);
r = roi - bg1;
if twoPass
    rv = r(valid);
    sg = 1.4826 * median(abs(rv - median(rv)));
    warm = r > 3 * sg;
    warm(~valid) = false;
    if grow > 0, warm = conv2(double(warm), ones(2*grow+1), 'same') > 0; end
    r2 = roi;  r2(warm) = NaN;
    [bg2, w2] = blurbg(r2, s);
    bg2(w2 < 0.02) = bg1(w2 < 0.02);        % too little plate left: keep the pass-1 background
    r = roi - bg2;
end
r(~valid) = NaN;
end

function [bg, w] = blurbg(roi, s)
valid = ~isnan(roi);
z = roi;  z(~valid) = 0;
w  = imgaussfilt(double(valid), s);
bg = imgaussfilt(z, s) ./ max(w, eps);
end

function roi = loadROI(folder, k)
fpath = fullfile(folder, sprintf('IV_%05d.xlsx', k));
if ~isfile(fpath), error('File not found: %s', fpath); end
[roi, ~] = extract_plate_roi(load_thermal_matrix(fpath), false, '', true);
end
