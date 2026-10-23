% CHECK_ROI_FIT  Checks that the fixed plate ROI (extract_plate_roi, outline A/B/C) sits on the plate in every frame of each folder.
% Flags frames with background or cold edges as CHECK, assigns each frame to a plate-temperature band and records the outline used.
% Edit the settings and run. Output: stats\frame_bands_<foldername>.csv, one per folder (read by corrosion_assess / painted_assess).

%Settings
baseDir = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\FINAL PIPELINES\paint pp';
%baseDir = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\painted pipeline';
addpath(baseDir);
folderNames = {'PAINT15','PAINT152','PAINT50','PAINT502','PAINT100','PAINT1002'};   % folders inside baseDir to check (one CSV each)
pattern     = 'IV_*.xlsx';     % files to check ('IV_013*.xlsx' for a subset)
stripPx     = 6;               % width of the edge strip checked, in pixels
coldDelta   = 4.5;             % minimum degC below the plate median that counts as "background"
coldK       = 2;               % ...or coldK * plate sigma below the median, whichever is larger
stripLimit  = 0.2;             % flag a side if more than this share of its strip is background
medDrop     = 3;               % flag a side if its strip median is this many degC below the ROI median (minimum limit)
medDropK    = 1.5;             % ...or medDropK * plate sigma, whichever is larger
coldLimit   = 1.0;             % flag the frame if more than this % of the ROI is cold
showFlagged = true;            % open the check figure for flagged frames
bandRanges  = [20 22; 22 26;   % [low high] degC, a frame is in a band if low <= plate median < high
               28 30; 30 32;
               38 40; 40 42;
               48 50; 50 52];
statsFolder = fullfile(baseDir, 'stats');

bandLabels = strings(size(bandRanges, 1), 1);
for j = 1:size(bandRanges, 1)
    bandLabels(j) = sprintf('%g-%g', bandRanges(j,1), bandRanges(j,2));
end
if ~exist(statsFolder, 'dir'), mkdir(statsFolder); end

for fi = 1:numel(folderNames)
    plateName = folderNames{fi};   % used in the CSV name
    folder    = fullfile(baseDir, plateName);
    files = dir(fullfile(folder, pattern));
    if isempty(files)
        warning('No files matching %s in %s - skipped.', pattern, folder);
        continue
    end

    %1.Check frames
    fprintf('\nChecking %d frames in %s\n', numel(files), folder);
    fprintf('%-16s %-8s %-8s %-8s %-8s %-26s %-10s %-8s %-10s %-4s %s\n', 'Frame', 'Median', 'Std', 'Cold%', 'Hot%', 'Worst edge strip (bg %)', 'Drop(degC)', 'Mean', 'Band', 'Pos', 'Status');

    nFlag = 0;
    nF = numel(files);
    bandOf = strings(nF, 1);  usable = false(nF, 1);
    rec = struct('frame', zeros(nF,1), 'median', zeros(nF,1), 'mean', zeros(nF,1), 'std', zeros(nF,1), ...
                 'edgeBgPct', zeros(nF,1), 'edgeDropC', zeros(nF,1), 'hotAreaPct', zeros(nF,1), 'status', strings(nF,1), 'pos', strings(nF,1));
    for i = 1:nF
        fpath = fullfile(folder, files(i).name);
        d = load_thermal_matrix(fpath);
        [~, info] = extract_plate_roi(d, false, files(i).name, true);   % quiet

        r = info.rows;  c = info.cols;  med = info.median;
        strips = {'top',    d(r(1):r(1)+stripPx-1,       c);
                  'bottom', d(r(end)-stripPx+1:r(end),   c);
                  'left',   d(r,                         c(1):c(1)+stripPx-1);
                  'right',  d(r,                         c(end)-stripPx+1:c(end))};
        cutoff = med - max(coldDelta, coldK * info.sigma);   % scale-aware "background" level (robust sigma)
        frac = zeros(1, 4);  drop = zeros(1, 4);
        for k = 1:4
            s = strips{k, 2};
            frac(k) = mean(s(:) < cutoff);
            drop(k) = med - median(s(:), 'omitnan');
        end
        [worst, wk] = max(frac);
        [worstDrop, dk] = max(drop);

        bad = {};
        if worst > stripLimit, bad{end+1} = sprintf('%s edge', strips{wk,1}); end %#ok<SAGROW>
        dropLimit = max(medDrop, medDropK * info.sigma);   % scale-aware limit on how much cooler an edge strip may be
        if worstDrop > dropLimit && ~(worst > stripLimit && wk == dk)
            bad{end+1} = sprintf('%s edge cooler', strips{dk,1}); %#ok<SAGROW>
        end
        if info.pctCold > coldLimit, bad{end+1} = 'many cold pixels'; end %#ok<SAGROW>
        if isempty(bad), status = 'OK'; else, status = ['CHECK: ' strjoin(bad, ', ')]; nFlag = nFlag + 1; end

        % band from the plate median
        bandStr = 'OUT';
        for j = 1:size(bandRanges, 1)
            lo = bandRanges(j,1);  hi = bandRanges(j,2);
            inBand = info.median >= lo && (info.median < hi || (j == size(bandRanges,1) && info.median <= hi));
            if inBand, bandStr = char(bandLabels(j)); break; end
        end
        bandOf(i) = bandStr;  usable(i) = isempty(bad);
        rec.frame(i) = str2double(regexp(files(i).name, '\d+', 'match', 'once'));
        rec.median(i) = info.median;  rec.mean(i) = info.mean;  rec.std(i) = info.std;
        rec.edgeBgPct(i) = 100*worst;  rec.edgeDropC(i) = worstDrop;
        rec.hotAreaPct(i) = info.hotAreaPct;  rec.status(i) = string(status);  rec.pos(i) = string(info.position);

        fprintf('%-16s %-8.2f %-8.2f %-8.2f %-8.2f %-26s %-10.2f %-8.2f %-10s %-4s %s\n', files(i).name, info.median, info.std, ...
            info.pctCold, info.hotAreaPct, sprintf('%s %.1f%%', strips{wk,1}, 100*worst), worstDrop, info.mean, bandStr, info.position, status);

        if showFlagged && ~isempty(bad)
            extract_plate_roi(d, true, files(i).name, true);
        end
    end

    %2.Summary
    fprintf('\n[%s] Frames per band (usable = status OK):\n', plateName);
    for lbl = [bandLabels; "OUT"]'
        inBand = (bandOf == lbl);
        names  = {files(usable & inBand).name};
        fprintf('  %-8s usable %2d / total %2d', lbl, numel(names), nnz(inBand));
        if ~isempty(names)
            fprintf(':  %s', strjoin(regexprep(names, '^IV_0*|\.xlsx$', ''), ' '));
        end
        fprintf('\n');
    end

    fprintf('Outline used:');
    for pn = ["A" "B" "C"], fprintf('  %s %d', pn, nnz(rec.pos == pn)); end
    fprintf('\n');

    %3.Save
    T = table(rec.frame, bandOf, rec.median, rec.mean, rec.std, rec.edgeBgPct, rec.edgeDropC, rec.hotAreaPct, rec.status, rec.pos, ...
        'VariableNames', {'frame','band','median','mean','std','edgeBgPct','edgeDropC','hotAreaPct','status','position'});
    T = sortrows(T, {'band','frame'});
    csvPath = fullfile(statsFolder, ['frame_bands_' plateName '.csv']);
    writetable(T, csvPath);
    fprintf('\nSaved band table: %s\n', csvPath);

    fprintf('%d of %d frames flagged in %s.\n', nFlag, numel(files), plateName);
end