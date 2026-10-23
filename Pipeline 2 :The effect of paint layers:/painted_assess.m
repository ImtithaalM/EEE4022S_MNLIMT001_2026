% PAINTED_ASSESS  Corrosion assessment of every coated plate in "plates", each against its own clean coated reference (same coat count).
% Frames flagged CHECK in stats\frame_bands_<folder>.csv (check_roi_fit) are excluded automatically. Run check_roi_fit first.
% Output: console tables, stats\assess_<plate>.csv (plus _frames.csv, _run.txt) and stats\painted_assess_summary.csv.

%Settings
%root        = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\painted pipeline';
root       = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\FINAL PIPELINES\paint pp';
addpath(root);
statsFolder = fullfile(root, 'stats');

plates = struct('name',     {'PAINT15','PAINT50','PAINT100','PAINT152','PAINT502','PAINT1002'}, ...
                'ref',      {'BARE PAINT','BARE PAINT','BARE PAINT','BARE PAINT2','BARE PAINT2','BARE PAINT2'}, ...
                'truthPct', {13.5, 55.3, 93, 13.5, 55.3, 93});

useRoiFlags     = true;     % exclude frames flagged CHECK in the check_roi_fit csv
skipFrames      = [];       % extra frame numbers never used (manual, on top of check_roi_fit flags)
printFrames     = false;    % true = also print the per-frame table (always saved to _frames.csv)
goodMinFrames   = 3;        % roll-up uses bands with at least this many reference and test frames
goodCloseness   = 3;        % roll-up uses bands with closeness at most this (degC)
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
frameSpreadMax  = 30;       % band is MIXED if per-frame areas span more than this many %-points
minUsableFrac   = 0.5;      % frame is usable if at least this share of ROI pixels are valid
bgBox     = [59 434 56 579];  % plate outer box [top bottom left right]; background is outside it
bgGap     = 15;               % px kept clear around the box
frameSize = [480 640];

if ~isfolder(statsFolder), mkdir(statsFolder); end
refCache = containers.Map('KeyType', 'char', 'ValueType', 'any');   % each reference is read once
allT = cell(1, numel(plates));

for tp = 1:numel(plates)
test = struct('name', plates(tp).name, 'folder', fullfile(root, plates(tp).name));
refName = plates(tp).ref;
refs = struct('name', {refName}, 'folder', {fullfile(root, refName)});
fprintf('\n################ TEST PLATE %s (reference %s) ################\n', test.name, refName);

%1.Inventory
if ~isKey(refCache, refName)
    R0 = inventory(refs.folder, [], inventoryStep, bgBox, bgGap, frameSize, minUsableFrac);
    fprintf('%s: %d frames read, %d usable.\n', refName, height(R0), nnz(R0.ok));
    R0 = addBand(R0, bandWidth, bandRange);
    if useRoiFlags, R0 = applyRoiFlags(R0, refs.folder, statsFolder, refName); end
    refCache(refName) = R0;
end
Inv.ref{1} = refCache(refName);
Inv.ref{1}.ok(ismember(Inv.ref{1}.frame, skipFrames)) = false;
Inv.test = inventory(test.folder, [], inventoryStep, bgBox, bgGap, frameSize, minUsableFrac);
fprintf('%s: %d frames read, %d usable.\n', test.name, height(Inv.test), nnz(Inv.test.ok));
Inv.test = addBand(Inv.test, bandWidth, bandRange);
if useRoiFlags, Inv.test = applyRoiFlags(Inv.test, test.folder, statsFolder, test.name); end
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
        if numel(iR) < minRef
            fprintf('band %s: skipped (%d reference frame(s), need %d; %d test frames).\n', b, numel(iR), minRef, numel(iT));
            continue
        end

        %2.Session match
        dBG = mean(Inv.test.bg(iT), 'omitnan') - mean(R.bg(iR), 'omitnan');
        dPB = mean(Inv.test.plateBG(iT), 'omitnan') - mean(R.plateBG(iR), 'omitnan');
        closeness = abs(dBG) + abs(dPB);

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
        if numel(unique(fv)) > 1, mixed = true;
        elseif max(fa) - min(fa) > frameSpreadMax, mixed = true; end
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
        rows(end+1, :) = {char(refs(r).name), char(b), numel(iR), numel(iT), dBG, dPB, closeness, lo, hi, sm, outside, ...
            char(verdict), sr, uaM, lg, char(branch), comb, char(cls)}; %#ok<SAGROW>
    end
end

if isempty(rows), warning('%s: no band had both test frames and at least %d reference frames. Skipped.', test.name, minRef); continue; end
T = cell2table(rows, 'VariableNames', {'Reference','Band','nRef','nTest','dBG','dPlateBG','Closeness','NullLo','NullHi', ...
    'StdRatio','OutsidePct','Verdict','SteelRefPct','UniformPct','LargestPct','Branch','AreaPct','Class'});
photoPct = plates(tp).truthPct;

%Results
fprintf('\n================ Corrosion detection assessment of "%s" ================\n', test.name);
hdr = sprintf('%-9s %-7s %-4s %-5s | %-9s | %-8s %-12s %-8s | %-8s %-8s | %s', ...
    'Ref', 'Band', 'nRef', 'nTest', 'Null', 'StdRatio', 'Verdict', 'Outside%', 'SteelRef', 'Area%', 'Class');
fprintf('%s\n%s\n', hdr, repmat('-', 1, numel(hdr)));
for k = 1:height(T)
    fprintf('%-9s %-7s %-4d %-5d | %-9s | %-8.2f %-12s %-8.0f | %-8.1f %-8.1f | %s\n', ...
        T.Reference{k}, T.Band{k}, T.nRef(k), T.nTest(k), sprintf('%.2f-%.2f', T.NullLo(k), T.NullHi(k)), ...
        T.StdRatio(k), T.Verdict{k}, T.OutsidePct(k), T.SteelRefPct(k), T.AreaPct(k), T.Class{k});
end
fprintf('%s\n', repmat('=', 1, numel(hdr)));
if ~isnan(photoPct), fprintf('Photo rust: %g %% of the plate\n', photoPct); end

TF = cell2table(frameRows, 'VariableNames', {'Reference','Band','Frame','StdRatio','Verdict','SteelRefPct','AreaPct','Branch'});
if printFrames
    fprintf('\n---- per frame ----\n');
    fprintf('%-9s %-7s %-7s | %-8s %-10s | %-8s %-8s %s\n', 'Ref', 'Band', 'Frame', 'StdRatio', 'Verdict', 'SteelRef', 'Area%', 'Branch');
    for k = 1:height(TF)
        fprintf('%-9s %-7s %-7d | %-8.2f %-10s | %-8.1f %-8.1f %s\n', TF.Reference{k}, TF.Band{k}, TF.Frame(k), TF.StdRatio(k), TF.Verdict{k}, TF.SteelRefPct(k), TF.AreaPct(k), TF.Branch{k});
    end
end

%Save
fid = fopen(fullfile(statsFolder, ['assess_' char(test.name) '_run.txt']), 'w');
fprintf(fid, 'painted_assess run %s\ntest: %s (%s)\nreference: %s (%s)\n', char(datetime('now')), char(test.name), test.folder, refName, refs.folder);
fprintf(fid, 'skipFrames: %s\nbandWidth %g, bandRange [%g %g], maxTestFrames %d, maxRefFrames %d, minRef %d, margin %g\n', mat2str(skipFrames), bandWidth, bandRange(1), bandRange(2), maxTestFrames, maxRefFrames, minRef, margin);
fprintf(fid, 'smoothSigma %g, pLow %g, pHigh %g, kSelf %g, gateHotPct %g, widespreadMin %g, localisedMin %g, frameSpreadMax %g\n\n', smoothSigma, pLow, pHigh, kSelf, gateHotPct, widespreadMin, localisedMin, frameSpreadMax);
fprintf(fid, '%s\n', runLog);  fclose(fid);
writetable(TF, fullfile(statsFolder, ['assess_' char(test.name) '_frames.csv']));
outFile = fullfile(statsFolder, ['assess_' char(test.name) '.csv']);
writetable(T, outFile);
fprintf('Saved %s\n', outFile);
allT{tp} = T;
end

%Roll-up (good bands only)
fprintf('\n================ Coated plates: roll-up ================\n');
hdr = sprintf('%-10s %-10s %-6s %-5s | %-18s | %-20s | %-24s | %s', 'Plate', 'Ref', 'Bands', 'Good', 'StdRatio med (rng)', 'Area% med (range)', 'Class (most common)', 'Photo%');
fprintf('%s\n%s\n', hdr, repmat('-', 1, numel(hdr)));
S = {};
for tp = 1:numel(plates)
    T = allT{tp};
    if isempty(T), continue; end
    good = T.nRef >= goodMinFrames & T.nTest >= goodMinFrames & T.Closeness <= goodCloseness;
    if ~any(good), good = true(height(T), 1); end
    cls = string(T.Class(good));
    [u, ~, ic] = unique(cls);  [~, im] = max(accumarray(ic, 1));  topClass = u(im);
    a = T.AreaPct(good);  s = T.StdRatio(good);
    fprintf('%-10s %-10s %-6d %-5d | %5.2f (%4.2f-%4.2f)   | %5.1f (%4.1f-%4.1f)      | %-24s | %g\n', plates(tp).name, plates(tp).ref, height(T), nnz(good), ...
        median(s), min(s), max(s), median(a, 'omitnan'), min(a), max(a), topClass, plates(tp).truthPct);
    S(end+1, :) = {plates(tp).name, plates(tp).ref, height(T), nnz(good), median(s), min(s), max(s), median(a, 'omitnan'), min(a), max(a), char(topClass), plates(tp).truthPct}; %#ok<SAGROW>
end
fprintf('%s\n', repmat('=', 1, numel(hdr)));
if ~isempty(S)
    TS = cell2table(S, 'VariableNames', {'Plate','Reference','Bands','GoodBands','StdMed','StdMin','StdMax','AreaMed','AreaMin','AreaMax','Class','PhotoPct'});
    writetable(TS, fullfile(statsFolder, 'painted_assess_summary.csv'));
    fprintf('Saved %s\n', fullfile(statsFolder, 'painted_assess_summary.csv'));
end

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
