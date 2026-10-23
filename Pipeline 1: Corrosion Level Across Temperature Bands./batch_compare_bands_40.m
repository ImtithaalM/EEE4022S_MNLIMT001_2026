% BATCH_COMPARE_BANDS  Compares test-plate frames (FULL anchor, light, HALF) against bare frames of the SAME temperature band and baseline source.
% Baselines are never pooled; all plates use bare2 = BARE NEW (2115-2123). Reads stats\frame_bands_<folder>.csv (run check_roi_fit first).
% THIS VERSION: bands 38-42 (38-40 and 40-42). Output: stats\band_results_38-42.csv

% ---- data folders ----
baseRoot   = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\matlab pipeline';
%baseRoot    = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\matlab pipeline';
addpath(baseRoot);
bareFolder  = fullfile(baseRoot, 'bare');
corrFolder  = fullfile(baseRoot, 'FULL');   % not '100' (painted)
statsFolder = fullfile(baseRoot, 'stats');
bandSet     = {'38-40', '40-42'};       % bands compared
bandTag     = '38-42';                  % output file name tag
if ~isfolder(statsFolder), mkdir(statsFolder); end

%1.Baselines
% Old bare (1337-1339) is a warmer-room session (BG 24.89) and is not used; 'BARE partial' is never a baseline.
bareSources = struct( ...
    'tag',    {'bare', 'bare2'}, ...
    'folder', {fullfile(baseRoot, 'bare'), fullfile(baseRoot, 'BARE NEW')}, ...
    'csv',    {fullfile(statsFolder, 'frame_bands_bare.csv'), fullfile(statsFolder, 'frame_bands_BARE NEW.csv')}, ...
    'keep',   {1337:1339, 2115:2123});
useSources = {'bare2'};     % add 'bare' only as a cross-check
minBase    = 2;             % skip a band if its baseline has fewer usable frames
bareSources = bareSources(ismember({bareSources.tag}, useSources));

comps = struct('name', {}, 'base', {}, 'null', {}, 'target', {}, 'loo', {}, 'light', {}, ...
               'tfolder', {}, 'trole', {}, 'bfolder', {});

%2.Test plates
% Only usable (status OK) frames. FULL ('100% corroded') sets the uniformity-index anchor.
% 3+ baseline frames: the last is held out as a bare null; 2+: a leave-one-out check is added.
testPlates = struct('name',      {'FULL', 'light', 'HALF'}, ...
                    'folder',    {fullfile(baseRoot, 'FULL'), fullfile(baseRoot, 'light'), fullfile(baseRoot, 'HALF')}, ...
                    'role',      {'100% corroded', 'LIGHT', 'HALF'}, ...
                    'baselines', {{'bare2'}, {'bare2'}, {'bare2'}});

%3.Session table
%   name        BG mean  Plate-BG
sessionTbl = {
    'FULL'   , 21.10, 19.28;     % main set 1746-1751
    'light'  , 20.24, 19.66;     % 1802-1813
    'HALF'   , 20.38, 19.52;     % main set 1610-1631
    'bare'   , 24.89, 14.55;     % 1337-1339
    'bare2'  , 20.36, 19.24};    % BARE NEW 2115-2123
fprintf('Session closeness (|dBG| + |dPlate-BG|, smaller = better baseline):\n');
for i = 1:3
    di = sessionTbl{i,2}; pi_ = sessionTbl{i,3};
    fprintf('  %-6s  vs bare %.2f   vs bare2 %.2f\n', sessionTbl{i,1}, ...
        abs(di - sessionTbl{4,2}) + abs(pi_ - sessionTbl{4,3}), abs(di - sessionTbl{5,2}) + abs(pi_ - sessionTbl{5,3}));
end

%4.Build comparisons
Tt_all  = cell(1, numel(testPlates));
okT_all = cell(1, numel(testPlates));
for p = 1:numel(testPlates)
    csvT = fullfile(statsFolder, ['frame_bands_' testPlates(p).name '.csv']);
    if ~isfile(csvT)
        if testPlates(p).name == "HALF"
            fprintf('HALF skipped: %s not found (run check_roi_fit on the HALF folder).\n', csvT);  continue
        end
        error('Band CSV missing: %s (run check_roi_fit on the %s folder).', csvT, testPlates(p).name);
    end
    Tt_all{p}  = readtable(csvT, 'TextType', 'string');
    okT_all{p} = startsWith(string(Tt_all{p}.status), "OK");
end

looDone = strings(0, 1);
for s = 1:numel(bareSources)
    src = bareSources(s);
    if ~isfile(src.csv), error('Band CSV missing: %s (run check_roi_fit on the %s folder).', src.csv, src.tag); end
    Tb  = readtable(src.csv, 'TextType', 'string');
    okB = startsWith(string(Tb.status), "OK");
    if ~isempty(src.keep), okB = okB & ismember(Tb.frame, src.keep); end
    for p = 1:numel(testPlates)
        tp = testPlates(p);  Tt = Tt_all{p};  okT = okT_all{p};
        if isempty(Tt), continue; end
        if ~ismember(src.tag, tp.baselines), continue; end
        bands = unique(strtrim(Tt.band(okT)), 'stable');
        for bi = 1:numel(bands)
            bnd = bands(bi);
            if contains(bnd, 'OUT') || ~any(bnd == string(bandSet)), continue; end
            baseFr = sort(Tb.frame(strtrim(Tb.band) == bnd & okB))';
            testFr = sort(Tt.frame(strtrim(Tt.band) == bnd & okT))';
            if numel(baseFr) < minBase || isempty(testFr)
                fprintf('%s vs %s band %s skipped: %d baseline (need %d) / %d %s usable frames.\n', ...
                    tp.name, src.tag, bnd, numel(baseFr), minBase, numel(testFr), tp.name);
                continue
            end
            allBare = baseFr;
            nullFr = [];
            if numel(baseFr) >= 3, nullFr = baseFr(end); baseFr = baseFr(1:end-1); end
            comps(end+1) = struct('name', sprintf('%s vs %s %s', tp.name, src.tag, bnd), ...
                'base', baseFr, 'null', nullFr, 'target', [], 'loo', false, 'light', testFr, ...
                'tfolder', tp.folder, 'trole', tp.role, 'bfolder', src.folder); %#ok<SAGROW>
            key = src.tag + "|" + bnd;
            if numel(allBare) >= 2 && ~ismember(key, looDone)    % one bare-vs-bare check per baseline and band
                comps(end+1) = struct('name', sprintf('band LOO %s %s', src.tag, bnd), ...
                    'base', allBare, 'null', [], 'target', [], 'loo', true, 'light', [], ...
                    'tfolder', '', 'trole', '', 'bfolder', src.folder); %#ok<SAGROW>
                looDone(end+1, 1) = key; %#ok<SAGROW>
            end
        end
    end
end

%5.Run comparisons
fields = {'offset','dT_mean','dT_max','dT_min','T_tilde', ...
          'N_tilde','N_tilde_med','std_ratio','iqr_ratio', ...
          'std_ratio_dt','iqr_ratio_dt'};
hdr = {'Comparison','Role','Frame','nBase','BaseTemp','FrameTemp','hotAreaPct'};
for f = 1:numel(fields)
    hdr{end+1} = [fields{f} '_mean']; %#ok<SAGROW>
    hdr{end+1} = [fields{f} '_sd'];   %#ok<SAGROW>
end
rows = {};

for c = 1:numel(comps)
    cp = comps(c);
    fprintf('Comparison %d/%d: %s\n', c, numel(comps), cp.name);
    if numel(cp.base) < 2
        fprintf('  note: only %d baseline frame, so the sd columns are not meaningful.\n', numel(cp.base));
    end

    baseROIs = cell(numel(cp.base), 1);
    meds = zeros(numel(cp.base), 1);
    for i = 1:numel(cp.base)
        baseROIs{i} = loadROI(cp.bfolder, cp.base(i));
        meds(i) = median(baseROIs{i}(:), 'omitnan');
    end
    baseTemp = mean(meds);

    if cp.loo
        % leave-one-out nulls
        for j = 1:numel(cp.base)
            others = setdiff(1:numel(cp.base), j);
            out = compute_indicators_multi(baseROIs(others), baseROIs{j}, '', true);
            [~, hj] = loadROI(cp.bfolder, cp.base(j));
            row = {cp.name, 'null (LOO)', cp.base(j), numel(others), mean(meds(others)), meds(j), hj};
            for f = 1:numel(fields)
                row{end+1} = out.mean.(fields{f}); %#ok<SAGROW>
                row{end+1} = out.sd.(fields{f});   %#ok<SAGROW>
            end
            rows(end+1, :) = row; %#ok<SAGROW>
        end
        continue
    end

    % held-out bare frames (null), then test-plate frames
    items = {};
    for k = cp.null,   items(end+1,:) = {'null (bare)', k, cp.bfolder}; end %#ok<SAGROW>
    for k = cp.light,  items(end+1,:) = {cp.trole, k, cp.tfolder}; end %#ok<SAGROW>

    for it = 1:size(items, 1)
        [roi, hA] = loadROI(items{it,3}, items{it,2});
        out = compute_indicators_multi(baseROIs, roi, '', true);
        row = {cp.name, items{it,1}, items{it,2}, numel(cp.base), baseTemp, median(roi(:), 'omitnan'), hA};
        for f = 1:numel(fields)
            row{end+1} = out.mean.(fields{f}); %#ok<SAGROW>
            row{end+1} = out.sd.(fields{f});   %#ok<SAGROW>
        end
        rows(end+1, :) = row; %#ok<SAGROW>
    end
end

if isempty(rows), error('No comparisons produced results. Check the band CSVs, bareSources and minBase.'); end
T = cell2table(rows, 'VariableNames', hdr);

%6.Uniformity index
% 0 = like bare, 1 = like 100% corroded. Anchor = median std ratio of FULL scored against the SAME baseline and band (NaN if none).
isC = strcmp(T.Role, '100% corroded');
grp = regexprep(T.Comparison, '^(FULL|light|HALF) vs ', '');
grp = regexprep(grp, '^band LOO ', '');
T.Baseline = grp;
T.uniformity_index = nan(height(T), 1);
T.uniformity_index_dt = nan(height(T), 1);      % plate slope removed
groups = unique(grp, 'stable');
for g = 1:numel(groups)
    m  = strcmp(grp, groups{g});
    mc = m & isC;
    if ~any(mc)
        fprintf('No FULL anchor for "%s": uniformity index left as NaN.\n', groups{g});
        continue
    end
    anchor   = median(T.std_ratio_mean(mc));
    anchorDT = median(T.std_ratio_dt_mean(mc));
    T.uniformity_index(m)    = uniformity_index(T.std_ratio_mean(m), anchor);
    T.uniformity_index_dt(m) = uniformity_index(T.std_ratio_dt_mean(m), anchorDT);
end
sel = {'Comparison','Baseline','Role','Frame','BaseTemp','FrameTemp','offset_mean','hotAreaPct', ...
       'std_ratio_mean','std_ratio_sd','uniformity_index','std_ratio_dt_mean','uniformity_index_dt','iqr_ratio_mean','dT_mean_mean','dT_mean_sd', ...
       'N_tilde_med_mean','N_tilde_med_sd'};

%7.Results
fprintf('\n%-38s %-14s %-6s %-9s %-10s %-8s %s\n', 'Comparison', 'Role', 'Frame', 'BaseTemp', 'FrameTemp', 'offset', '<strong>STD_RATIO</strong>');
for r = 1:height(T)
    fprintf('%-38s %-14s %-6d %-9.2f %-10.2f %-8.2f <strong>%.3f</strong>  (sd %.3f)\n', ...
        T.Comparison{r}, T.Role{r}, T.Frame(r), T.BaseTemp(r), T.FrameTemp(r), T.offset_mean(r), ...
        T.std_ratio_mean(r), T.std_ratio_sd(r));
end
fprintf('\n');
disp(T(:, sel));

% std ratio by comparison and role (compare only within the same baseline)
names = unique(T.Comparison, 'stable');
fprintf('\nSTD_RATIO summary by comparison and role:\n');
fprintf('%-38s %-14s %-4s %-8s %-8s %-8s %-8s\n', 'Comparison', 'Role', 'n', 'mean', 'sd', 'min', 'max');
for n = 1:numel(names)
    m = strcmp(T.Comparison, names{n});
    roles = unique(T.Role(m), 'stable');
    for ri = 1:numel(roles)
        v = T.std_ratio_mean(m & strcmp(T.Role, roles{ri}));
        fprintf('%-38s %-14s %-4d %-8.3f %-8.3f %-8.3f %-8.3f\n', names{n}, roles{ri}, numel(v), ...
            mean(v, 'omitnan'), std(v, 'omitnan'), min(v), max(v));
    end
end
fprintf('\n');

%8.Save
outName = ['band_results_' bandTag '.csv'];
writetable(T, fullfile(statsFolder, outName));
fprintf('Saved %s\n', fullfile(statsFolder, outName));

% ---------------------------------------------------------------------
function [roi, hotAreaPct] = loadROI(folder, k)
% Plate ROI of frame k (specks masked, large hot regions kept) and its hotAreaPct.
fpath = fullfile(folder, sprintf('IV_%05d.xlsx', k));
if ~isfile(fpath), error('File not found: %s', fpath); end
d = load_thermal_matrix(fpath);
[roi, info] = extract_plate_roi(d, false, '', true);
hotAreaPct = info.hotAreaPct;
end
