% CHECK_ROI_FIT  Checks the fixed plate ROI (extract_plate_roi) still sits on the plate for every frame in each listed folder.
% Flags a side if its edge strip is mostly background or much cooler than the plate, and assigns each frame to a temperature band.
% Edit the settings and run. Output: stats\frame_bands_<foldername>.csv, one per folder.

% ---- base path ----
baseDir = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\FINAL PIPELINES\corr pp';
%baseDir = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\matlab pipeline';
addpath(baseDir);

% ---- settings ----
folderNames = {'BARE NEW' ,'FULL'};
% folderNames = {'bare', 'FULL', 'light', 'BARE partial'};
pattern     = 'IV_*.xlsx';   % files to check ('IV_013*.xlsx' for a subset)
stripPx     = 6;             % width of the edge strip checked (px)
coldDelta   = 4.5;           % min degC below plate median that counts as background
coldK       = 2;             % or coldK * plate sigma below median, whichever is larger
stripLimit  = 0.2;           % flag a side if more than this fraction of its strip is background
medDrop     = 3;             % flag a side if its strip median is this many degC below the ROI median (minimum)
medDropK    = 1.5;           % or medDropK * plate sigma, whichever is larger
coldLimit   = 1.0;           % flag the frame if more than this % of the ROI is cold-flagged
showFlagged = true;          % open the check figure for flagged frames
% Temperature bands (degC), one [low high] row each: low <= plate median < high
bandRanges  = [20 22; 22 26;
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
    plateName = folderNames{fi};
    folder    = fullfile(baseDir, plateName);
    files = dir(fullfile(folder, pattern));
    if isempty(files)
        warning('No files matching %s in %s - skipped.', pattern, folder);
        continue
    end

    fprintf('\nChecking %d frames in %s\n', numel(files), folder);
    fprintf('%-16s %-8s %-8s %-8s %-8s %-26s %-10s %-8s %-10s %s\n', 'Frame', 'Median', 'Std', 'Cold%', 'Hot%', 'Worst edge strip (bg %)', 'Drop(degC)', 'Mean', 'Band', 'Status');

    nFlag = 0;
    nF = numel(files);
    bandOf = strings(nF, 1);  usable = false(nF, 1);
    rec = struct('frame', zeros(nF,1), 'median', zeros(nF,1), 'mean', zeros(nF,1), 'std', zeros(nF,1), ...
                 'edgeBgPct', zeros(nF,1), 'edgeDropC', zeros(nF,1), 'hotAreaPct', zeros(nF,1), 'status', strings(nF,1));
    for i = 1:nF
        fpath = fullfile(folder, files(i).name);
        d = load_thermal_matrix(fpath);
        [~, info] = extract_plate_roi(d, false, files(i).name, true);

        %1.Edge strips
        r = info.rows;  c = info.cols;  med = info.median;
        strips = {'top',    d(r(1):r(1)+stripPx-1,       c);
                  'bottom', d(r(end)-stripPx+1:r(end),   c);
                  'left',   d(r,                         c(1):c(1)+stripPx-1);
                  'right',  d(r,                         c(end)-stripPx+1:c(end))};
        cutoff = med - max(coldDelta, coldK * info.sigma);
        frac = zeros(1, 4);  drop = zeros(1, 4);
        for k = 1:4
            s = strips{k, 2};
            frac(k) = mean(s(:) < cutoff);
            drop(k) = med - median(s(:), 'omitnan');
        end
        [worst, wk] = max(frac);
        [worstDrop, dk] = max(drop);

        %2.Flags
        bad = {};
        if worst > stripLimit, bad{end+1} = sprintf('%s edge', strips{wk,1}); end %#ok<SAGROW>
        dropLimit = max(medDrop, medDropK * info.sigma);
        if worstDrop > dropLimit && ~(worst > stripLimit && wk == dk)
            bad{end+1} = sprintf('%s edge cooler', strips{dk,1}); %#ok<SAGROW>
        end
        if info.pctCold > coldLimit, bad{end+1} = 'many cold pixels'; end %#ok<SAGROW>
        if isempty(bad), status = 'OK'; else, status = ['CHECK: ' strjoin(bad, ', ')]; nFlag = nFlag + 1; end

        %3.Band
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
        rec.hotAreaPct(i) = info.hotAreaPct;  rec.status(i) = string(status);

        fprintf('%-16s %-8.2f %-8.2f %-8.2f %-8.2f %-26s %-10.2f %-8.2f %-10s %s\n', files(i).name, info.median, info.std, ...
            info.pctCold, info.hotAreaPct, sprintf('%s %.1f%%', strips{wk,1}, 100*worst), worstDrop, info.mean, bandStr, status);

        if showFlagged && ~isempty(bad)
            extract_plate_roi(d, true, files(i).name, true);
        end
    end

    %4.Summary
    fprintf('\n[%s] Frames per band (usable = status OK, total = all frames in the band):\n', plateName);
    for lbl = [bandLabels; "OUT"]'
        inBand = (bandOf == lbl);
        names  = {files(usable & inBand).name};
        fprintf('  %-8s usable %2d / total %2d', lbl, numel(names), nnz(inBand));
        if ~isempty(names)
            fprintf(':  %s', strjoin(regexprep(names, '^IV_0*|\.xlsx$', ''), ' '));
        end
        fprintf('\n');
    end

    %5.Save
    T = table(rec.frame, bandOf, rec.median, rec.mean, rec.std, rec.edgeBgPct, rec.edgeDropC, rec.hotAreaPct, rec.status, ...
        'VariableNames', {'frame','band','median','mean','std','edgeBgPct','edgeDropC','hotAreaPct','status'});
    T = sortrows(T, {'band','frame'});
    csvPath = fullfile(statsFolder, ['frame_bands_' plateName '.csv']);
    writetable(T, csvPath);
    fprintf('\nSaved band table: %s\n', csvPath);
    fprintf('%d of %d frames flagged in %s.\n', nFlag, numel(files), plateName);
end
