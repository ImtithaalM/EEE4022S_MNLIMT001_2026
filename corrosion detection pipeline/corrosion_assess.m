% CORROSION_ASSESS  Assesses one folder of thermal frames against one or two reference folders, band by band (2 C plate-temperature bands).
% Frames flagged CHECK in stats\frame_bands_<folder>.csv (check_roi_fit) are excluded automatically. Run check_roi_fit first.
% Output: console tables, a std-ratio vs temperature plot, and stats\assess_<testName>.csv (plus _frames.csv, _run.txt, .png).

% ---- settings ----
root        = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\FINAL PIPELINES\corr pp';
%root       = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\matlab pipeline';
addpath(root);
statsFolder = fullfile(root, 'stats');

test = struct('name', 'FULL', 'folder', fullfile(root, 'FULL/'), 'frames', []);   % frames [] = every frame
refs = struct('name',   {'bare2'}, ...
              'folder', {fullfile(root, 'BARE NEW')}, ...
              'frames', {[]}, ...
              'exclude', {[2146 2152 2153 2123]});   % exclude = frames not used as reference; add a second entry for a second reference

bandWidth       = 2;        % degC per band
bandRange       = [20 60];  % plate-median range considered (degC)
inventoryStep   = 1;        % read every n-th frame in the inventory
maxTestFrames   = 12;       % test frames scored per band
maxRefFrames    = 12;       % reference frames used per band
minRef          = 2;        % band skipped if fewer usable reference frames
margin          = 0.05;     % std-ratio margin around the reference null range
smoothSigma     = 2;        % px smoothing before the steel-referenced area
pLow = 20;  pHigh = 95;     % steel level percentile, and percentile on the reference that sets DELTA
kSelf           = 3;        % uniform estimate: cooler than kSelf robust sigmas = steel-like
gateHotPct      = 5;        % steel-referenced area must be below this % for the uniform estimate
widespreadMin   = 30;       % PATCHY with area above this % = PARTIAL / WIDESPREAD (provisional)
localisedMin    = 7;        % PATCHY with area above this % = LOCALISED (provisional)
closenessWarn   = 3;        % degC: warn if |dBG| + |dPlate-BG| exceeds this (provisional)
skipFrames      = [];       % extra frame numbers never used (manual, on top of check_roi_fit flags)
useRoiFlags     = true;     % exclude frames flagged CHECK in the check_roi_fit csv
frameSpreadWarn = 30;       % warn if per-frame areas in a band span more than this many %-points
minUsableFrac   = 0.5;      % frame is usable if at least this share of ROI pixels are valid
bgBox     = [59 434 56 579];  % plate outer box [top bottom left right]; background is outside it
bgGap     = 15;               % px kept clear around the box
frameSize = [480 640];

if ~isfolder(statsFolder), mkdir(statsFolder); end

%1.Inventory
fprintf('\n=== 1. FRAME INVENTORY ===\n');
Inv.test = inventory(test.folder, test.frames, inventoryStep, bgBox, bgGap, frameSize, minUsableFrac);
fprintf('%s: %d frames read, %d usable.\n', test.name, height(Inv.test), nnz(Inv.test.ok));
for r = 1:numel(refs)
    Inv.ref{r} = inventory(refs(r).folder, refs(r).frames, inventoryStep, bgBox, bgGap, frameSize, minUsableFrac);
    fprintf('%s: %d frames read, %d usable.\n', refs(r).name, height(Inv.ref{r}), nnz(Inv.ref{r}.ok));
end
Inv.test = addBand(Inv.test, bandWidth, bandRange);
if useRoiFlags, Inv.test = applyRoiFlags(Inv.test, test.folder, statsFolder, test.name); end
for r = 1:numel(refs)
    Inv.ref{r} = addBand(Inv.ref{r}, bandWidth, bandRange);
    if useRoiFlags, Inv.ref{r} = applyRoiFlags(Inv.ref{r}, refs(r).folder, statsFolder, refs(r).name); end
    Inv.ref{r}.ok(ismember(Inv.ref{r}.frame, refs(r).exclude)) = false;
    Inv.ref{r}.ok(ismember(Inv.ref{r}.frame, skipFrames)) = false;
end
Inv.test.ok(ismember(Inv.test.frame, skipFrames)) = false;

bandsTest = unique(Inv.test.band(Inv.test.ok & Inv.test.band ~= ""), 'stable');
fprintf('\nTest plate frames per band (usable): ');
for b = bandsTest', fprintf('%s:%d  ', b, nnz(Inv.test.ok & Inv.test.band == b)); end
fprintf('\n');

rows = {};  frameRows = {};  runLog = strings(0, 1);
for r = 1:numel(refs)
    R = Inv.ref{r};
    fprintf('\n############ reference "%s" ############\n', refs(r).name);
    for b = sort(bandsTest)'
        iT = find(Inv.test.ok & Inv.test.band == b);
        iR = find(R.ok & R.band == b);
        warn = strings(0, 1);
        if numel(iR) < minRef
            fprintf('band %s: skipped (%d reference frame(s), need %d; %d test frames).\n', b, numel(iR), minRef, numel(iT));
            continue
        end
        if numel(iR) < 4, warn(end+1) = "thin baseline (" + numel(iR) + " ref frames)"; end %#ok<SAGROW>
        if numel(iT) < 3, warn(end+1) = "few test frames (" + numel(iT) + ")"; end %#ok<SAGROW>

        %2.Session match
        dBG = mean(Inv.test.bg(iT), 'omitnan') - mean(R.bg(iR), 'omitnan');
        dPB = mean(Inv.test.plateBG(iT), 'omitnan') - mean(R.plateBG(iR), 'omitnan');
        closeness = abs(dBG) + abs(dPB);
        if closeness > closenessWarn, warn(end+1) = sprintf("session mismatch (closeness %.1f)", closeness); end %#ok<SAGROW>

        % frames actually scored
        iT = pick(iT, maxTestFrames);  iR = pick(iR, maxRefFrames);
        refROI = cell(numel(iR), 1);
        for i = 1:numel(iR), refROI{i} = loadROI(refs(r).folder, R.frame(iR(i))); end
        tROI = cell(numel(iT), 1);
        for i = 1:numel(iT), tROI{i} = loadROI(test.folder, Inv.test.frame(iT(i))); end

        %3.Detection
        nullR = nan(numel(iR), 1);
        for j = 1:numel(iR)
            others = setdiff(1:numel(iR), j);
            o = compute_indicators_multi(refROI(others), refROI{j}, '', true);
            nullR(j) = o.mean.std_ratio;
        end
        stdT = nan(numel(iT), 1);
        for i = 1:numel(iT)
            o = compute_indicators_multi(refROI, tROI{i}, '', true);
            stdT(i) = o.mean.std_ratio;
        end
        lo = min(nullR);  hi = max(nullR);
        sm = mean(stdT, 'omitnan');
        if sm < lo - margin, verdict = "UNIFORM"; elseif sm > hi + margin, verdict = "PATCHY"; else, verdict = "BARE-LIKE"; end
        outside = 100 * mean(stdT < lo - margin | stdT > hi + margin);

        %4.Sizing
        dl = nan(numel(iR), 1);
        for j = 1:numel(iR)
            v = smoothVec(refROI{j}, smoothSigma);
            dl(j) = prctile(v, pHigh) - prctile(v, pLow);
        end
        DELTA = median(dl);
        e3a = nan(numel(iT), 1);  ua = e3a;  largest = e3a;
        for i = 1:numel(iT)
            [v, tS, mask] = smoothVec(tROI{i}, smoothSigma);
            p20 = prctile(v, pLow);
            e3a(i) = 100 * mean(v > p20 + DELTA);
            hot = mask & (tS > p20 + DELTA);
            cc = bwconncomp(hot);
            if cc.NumObjects > 0, largest(i) = 100 * max(cellfun(@numel, cc.PixelIdxList)) / nnz(mask); else, largest(i) = 0; end
            ua(i) = 100 - coldShare(tS, mask, kSelf);
        end
        % per-frame verdict and area
        fv = strings(numel(iT), 1);  fa = nan(numel(iT), 1);  fbr = strings(numel(iT), 1);
        for i = 1:numel(iT)
            if stdT(i) < lo - margin, fv(i) = "UNIFORM"; elseif stdT(i) > hi + margin, fv(i) = "PATCHY"; else, fv(i) = "BARE-LIKE"; end
            if fv(i) == "UNIFORM" && e3a(i) < gateHotPct, fa(i) = ua(i); fbr(i) = "uniform"; else, fa(i) = e3a(i); fbr(i) = "contrast"; end
            if fv(i) == "BARE-LIKE", fa(i) = 0; end
            frameRows(end+1, :) = {char(refs(r).name), char(b), Inv.test.frame(iT(i)), stdT(i), char(fv(i)), e3a(i), fa(i), char(fbr(i))}; %#ok<SAGROW>
        end
        runLog(end+1) = sprintf("ref %s band %s | reference frames: %s | test frames: %s", refs(r).name, b, mat2str(R.frame(iR)'), mat2str(Inv.test.frame(iT)')); %#ok<SAGROW>
        mixed = false;
        if numel(unique(fv)) > 1, mixed = true; warn(end+1) = "frames disagree (" + strjoin(unique(fv)', "/") + ")"; %#ok<SAGROW>
        elseif max(fa) - min(fa) > frameSpreadWarn, mixed = true; warn(end+1) = sprintf("frame areas span %.0f-%.0f %%", min(fa), max(fa)); end %#ok<SAGROW>
        sr = mean(e3a, 'omitnan');  uaM = mean(ua, 'omitnan');  lg = mean(largest, 'omitnan');
        gate = (verdict == "UNIFORM") && sr < gateHotPct;
        if gate, comb = uaM; branch = "uniform"; else, comb = sr; branch = "contrast"; end

        %5.Class
        if verdict == "BARE-LIKE"
            cls = "BARE-LIKE";
        elseif verdict == "UNIFORM"
            if gate, cls = "UNIFORM"; else, cls = "UNIFORM (hot area present)"; end
        else
            if sr > widespreadMin, cls = "PARTIAL / WIDESPREAD";
            elseif sr > localisedMin, cls = "LOCALISED";
            else, cls = "UNEVEN"; end
        end
        if mixed, cls = "MIXED";  comb = NaN;  sr = NaN;  uaM = NaN;  lg = NaN; end
        if numel(warn) == 0, wtxt = ""; else, wtxt = strjoin(warn, "; "); end
        rows(end+1, :) = {char(refs(r).name), char(b), numel(iR), numel(iT), dBG, dPB, closeness, lo, hi, sm, outside, ...
            char(verdict), sr, uaM, lg, char(branch), comb, char(cls), char(wtxt)}; %#ok<SAGROW>
    end
end

if isempty(rows), error('No band had both test frames and at least %d reference frames. Check the folders, bandRange and minRef.', minRef); end
T = cell2table(rows, 'VariableNames', {'Reference','Band','nRef','nTest','dBG','dPlateBG','Closeness','NullLo','NullHi', ...
    'StdRatio','OutsidePct','Verdict','SteelRefPct','UniformPct','LargestPct','Branch','AreaPct','Class','Warnings'});

%Results
fprintf('\n================ ASSESSMENT OF "%s" ================\n', test.name);
fprintf('%-8s %-7s %-4s %-5s | %-9s %-6s | %-8s %-12s %-8s | %-8s %-8s %-9s | %s\n', ...
    'Ref', 'Band', 'nRef', 'nTest', 'Closeness', 'Null', 'StdRatio', 'Verdict', 'Outside%', 'SteelRef', 'Area%', 'Branch', 'Class');
for k = 1:height(T)
    fprintf('%-8s %-7s %-4d %-5d | %-9.2f %-6s | %-8.2f %-12s %-8.0f | %-8.1f %-8.1f %-9s | %s\n', ...
        T.Reference{k}, T.Band{k}, T.nRef(k), T.nTest(k), T.Closeness(k), sprintf('%.2f-%.2f', T.NullLo(k), T.NullHi(k)), ...
        T.StdRatio(k), T.Verdict{k}, T.OutsidePct(k), T.SteelRefPct(k), T.AreaPct(k), T.Branch{k}, T.Class{k});
    if ~isempty(T.Warnings{k}), fprintf('         ! %s\n', T.Warnings{k}); end
end
fprintf('======================================================================\n');

%Summary
for r = 1:numel(refs)
    m = strcmp(T.Reference, char(refs(r).name));
    if ~any(m), continue; end
    V = T.Verdict(m);
    fprintf('\n%s vs %s: %d bands | PATCHY %d | UNIFORM %d | BARE-LIKE %d | MIXED %d\n', test.name, refs(r).name, nnz(m), ...
        nnz(strcmp(V, 'PATCHY')), nnz(strcmp(V, 'UNIFORM')), nnz(strcmp(V, 'BARE-LIKE')), nnz(strcmp(T.Class(m), 'MIXED')));
    [x, o] = sort(str2double(extractBefore(string(T.Band(m)), '-')) + bandWidth/2);
    nl = T.NullLo(m) - margin;  nl = nl(o);
    nh = T.NullHi(m) + margin;  nh = nh(o);
    y  = T.StdRatio(m);         y  = y(o);
    figure('Name', sprintf('%s vs %s', test.name, refs(r).name));  hold on;
    fill([x; flipud(x)], [nl; flipud(nh)], [0.8 0.8 0.8], 'EdgeColor', 'none');
    plot(x, y, 'o-', 'LineWidth', 1.5);
    xlabel('Plate temperature (degC)');  ylabel('Std ratio');  grid on;
    title(sprintf('%s vs %s', test.name, refs(r).name));
    legend({'Reference null range', char(test.name)}, 'Location', 'best');
    saveas(gcf, fullfile(statsFolder, ['assess_' char(test.name) '_' char(refs(r).name) '.png']));
end

TF = cell2table(frameRows, 'VariableNames', {'Reference','Band','Frame','StdRatio','Verdict','SteelRefPct','AreaPct','Branch'});
fprintf('\n---- per frame ----\n');
fprintf('%-8s %-7s %-7s | %-8s %-10s | %-8s %-8s %s\n', 'Ref', 'Band', 'Frame', 'StdRatio', 'Verdict', 'SteelRef', 'Area%', 'Branch');
for k = 1:height(TF)
    fprintf('%-8s %-7s %-7d | %-8.2f %-10s | %-8.1f %-8.1f %s\n', TF.Reference{k}, TF.Band{k}, TF.Frame(k), TF.StdRatio(k), TF.Verdict{k}, TF.SteelRefPct(k), TF.AreaPct(k), TF.Branch{k});
end

%Save
fid = fopen(fullfile(statsFolder, ['assess_' char(test.name) '_run.txt']), 'w');
fprintf(fid, 'corrosion_assess run %s\ntest: %s (%s)\n', char(datetime('now')), char(test.name), test.folder);
for r = 1:numel(refs), fprintf(fid, 'reference: %s (%s) exclude: %s\n', char(refs(r).name), refs(r).folder, mat2str(refs(r).exclude)); end
fprintf(fid, 'skipFrames: %s\nbandWidth %g, bandRange [%g %g], maxTestFrames %d, maxRefFrames %d, minRef %d, margin %g\n', mat2str(skipFrames), bandWidth, bandRange(1), bandRange(2), maxTestFrames, maxRefFrames, minRef, margin);
fprintf(fid, 'smoothSigma %g, pLow %g, pHigh %g, kSelf %g, gateHotPct %g, widespreadMin %g, localisedMin %g, closenessWarn %g, frameSpreadWarn %g\n\n', smoothSigma, pLow, pHigh, kSelf, gateHotPct, widespreadMin, localisedMin, closenessWarn, frameSpreadWarn);
fprintf(fid, '%s\n', runLog);  fclose(fid);
writetable(TF, fullfile(statsFolder, ['assess_' char(test.name) '_frames.csv']));
outFile = fullfile(statsFolder, ['assess_' char(test.name) '.csv']);
writetable(T, outFile);
fprintf('Saved %s\n', outFile);

% ---------------------------------------------------------------------
function Inv = inventory(folder, wantFrames, step, bgBox, bgGap, frameSize, minUsableFrac)
% One row per frame: number, plate median, usable flag, hotAreaPct, background mean, plate minus background.
files = dir(fullfile(folder, 'IV_*.xlsx'));
if isempty(files), error('No IV_*.xlsx frames found in %s.', folder); end
fr = zeros(numel(files), 1);
for i = 1:numel(files)
    t = regexp(files(i).name, '^IV_(\d+)\.xlsx$', 'tokens', 'once');
    if ~isempty(t), fr(i) = str2double(t{1}); end
end
fr = sort(fr(fr > 0));
if ~isempty(wantFrames), fr = fr(ismember(fr, wantFrames)); end
fr = fr(1:step:end);
n = numel(fr);
med = nan(n, 1);  ok = false(n, 1);  hotA = nan(n, 1);  bg = nan(n, 1);
bgRows = 100:400;  bgLeft = 5:40;  bgRight = 600:635;   % background strips either side of the plate (as in background_check.m)
for i = 1:n
    fpath = fullfile(folder, sprintf('IV_%05d.xlsx', fr(i)));
    try
        d = load_thermal_matrix(fpath);
        [roi, info] = extract_plate_roi(d, false, '', true);
        v = ~isnan(roi);
        ok(i) = nnz(v) / numel(roi) >= minUsableFrac;
        med(i) = median(roi(v));
        hotA(i) = info.hotAreaPct;
        if isequal(size(d), frameSize)
            bg(i) = (median(reshape(d(bgRows, bgLeft), [], 1), 'omitnan') + median(reshape(d(bgRows, bgRight), [], 1), 'omitnan')) / 2;
        end
    catch ME
        fprintf('  frame %d skipped (%s)\n', fr(i), ME.message);
    end
    if mod(i, 25) == 0, fprintf('  %s: %d/%d frames\n', folder, i, n); end
end
Inv = table(fr, med, ok, hotA, bg, med - bg, 'VariableNames', {'frame', 'median', 'ok', 'hotAreaPct', 'bg', 'plateBG'});
end

function Inv = addBand(Inv, w, rng)
% Band label 'lo-hi' from the plate median; "" outside rng or without a median.
lo = floor(Inv.median / w) * w;
lab = compose('%d-%d', lo, lo + w);
lab = string(lab);
lab(isnan(Inv.median) | Inv.median < rng(1) | Inv.median >= rng(2)) = "";
Inv.band = lab;
end

function idx = pick(idx, nMax)
% Evenly spaced subset of at most nMax entries.
if numel(idx) > nMax, idx = idx(unique(round(linspace(1, numel(idx), nMax)))); end
end

function roi = loadROI(folder, k)
d = load_thermal_matrix(fullfile(folder, sprintf('IV_%05d.xlsx', k)));
[roi, ~] = extract_plate_roi(d, false, '', true);
end

function [v, tS, mask] = smoothVec(roi, sigma)
% Smoothed temperature: v = valid-pixel vector, tS = smoothed image (NaN outside plate), mask = valid pixels.
mask = ~isnan(roi);
z = roi;  z(~mask) = median(roi(mask));
tS = imgaussfilt(z, sigma);  tS(~mask) = NaN;
v = tS(mask);
end

function c = coldShare(tS, mask, kSelf)
% Percent of valid pixels cooler than kSelf robust sigmas of the plate's own spread (plane and median removed).
[X, Y] = meshgrid(1:size(tS, 2), 1:size(tS, 1));
t1 = tS(mask);
Xd = [X(mask) Y(mask) ones(nnz(mask), 1)];
r = t1 - Xd * (Xd \ t1);
r = r - median(r);
sg = 1.4826 * median(abs(r));
c = 100 * mean(r < -kSelf * sg);
end

function T = applyRoiFlags(T, folder, statsFolder, name)
% Sets ok = false for frames whose status in stats\frame_bands_<folder>.csv is not OK.
[~, fname] = fileparts(char(regexprep(string(folder), '[\\/]+$', '')));
csv = fullfile(statsFolder, ['frame_bands_' fname '.csv']);
if ~isfile(csv)
    fprintf('%s: %s not found, ROI flags not applied (run check_roi_fit).\n', name, csv);
    return
end
C = readtable(csv, 'TextType', 'string');
bad = C.frame(~startsWith(strtrim(C.status), "OK"));
hit = ismember(T.frame, bad) & T.ok;
T.ok(ismember(T.frame, bad)) = false;
fprintf('%s: %d frame(s) excluded by check_roi_fit flags.\n', name, nnz(hit));
end