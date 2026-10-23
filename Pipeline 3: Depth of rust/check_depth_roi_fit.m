% CHECK_DEPTH_ROI_FIT  Checks that the plate ROI and depth rectangles (extract_depth_roi / depth_geometry) still sit on the plate, and tabulates each active rectangle per band.
% Frames are banded by the reference median (plate outside the rectangles); dT = rectangle median - reference median. Flagged frames are CHECK, frames in no band are OUT.
% Edit the settings and run. Output: stats\depth_frame_bands_<foldername>.csv, one per folder (read by analyse_depth).

%Settings
%baseDir = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\rust pipeline';
baseDir = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\rust pipeline';
addpath(baseDir);
folderNames = {'2mm new', '4mm new'};   % folders inside baseDir to check (one CSV each)
pattern     = 'IV_*.xlsx';     % files to check ('IV_023*.xlsx' for a subset)
stripPx     = 6;               % width of the edge strip checked, in pixels
coldDelta   = 4.5;             % minimum degC below the reference median that counts as background
coldK       = 2;               % ...or coldK * reference sigma, whichever is larger
stripLimit  = 0.35;            % flag a side if more than this share of its strip is background (0.2 in check_roi_fit; the cold patch fills ~22% of the bottom strip)
medDrop     = 3;               % flag a side if its strip median is this many degC below the reference median (minimum limit)
medDropK    = 1.5;             % ...or medDropK * reference sigma, whichever is larger
stripSkipGrow = 60;            % px: edge-strip pixels near the inactive rectangles / refExclude patch are ignored (0 = off)
rectMinValid = 80;             % flag the frame if an active rectangle has under this % valid pixels
showFlagged = true;            % open the check figure for flagged frames
showPerBand = true;            % one check figure per band (first frame in the band)
bandRanges  = [22 26;   % [low high] degC, a frame is in a band if low <= reference median < high
               28 32;
               38 42;
               48 52];
statsFolder = fullfile(baseDir, 'stats');

bandLabels = strings(size(bandRanges, 1), 1);
for j = 1:size(bandRanges, 1)
    bandLabels(j) = sprintf('%g-%g', bandRanges(j,1), bandRanges(j,2));
end
if ~exist(statsFolder, 'dir'), mkdir(statsFolder); end

stripNames = {'top', 'right', 'bottom', 'left'};
fitPrinted = false;

for fi = 1:numel(folderNames)
    plateName = folderNames{fi};                  % used in the CSV name
    folder    = fullfile(baseDir, plateName);
    plateKey  = '2mm';                            % geometry chosen from the folder name
    if contains(lower(plateName), '4mm'), plateKey = '4mm'; end

    %1.Geometry
    geom       = depth_geometry(plateKey);
    rectNames  = geom.rectNames;                  % all rectangles
    nR         = numel(rectNames);
    act        = geom.activeIdx;                  % rectangles that get tabulated

    fitPrinted = false;                           % rectangle clearance printed once per folder
    files = dir(fullfile(folder, pattern));
    if isempty(files)
        warning('No files matching %s in %s - skipped.', pattern, folder);
        continue
    end

    %2.Check frames
    fprintf('\nChecking %d frames in %s  (active rectangles: %s)\n', numel(files), folder, strjoin(rectNames(act), ' '));
    fprintf('%-16s %-7s %-8s', 'Frame', 'Ref', 'Band');
    for k = act, fprintf(' %-7s', ['dT_' rectNames{k}]); end
    fprintf(' %-26s %s\n', 'Worst edge strip (bg %)', 'Status');

    nFlag = 0;
    nF = numel(files);
    bandOf = strings(nF, 1);  usable = false(nF, 1);
    rec = struct('frame', zeros(nF,1), 'ref', zeros(nF,1), 'sigma', zeros(nF,1), ...
                 'edgeBgPct', zeros(nF,1), 'edgeDropC', zeros(nF,1), 'status', strings(nF,1), ...
                 'rMedian', nan(nF,nR), 'rMean', nan(nF,nR), 'rStd', nan(nF,nR), ...
                 'rDT', nan(nF,nR), 'rValid', nan(nF,nR));
    for i = 1:nF
        fpath = fullfile(folder, files(i).name);
        d = load_thermal_matrix(fpath);
        [~, info] = extract_depth_roi(d, false, files(i).name, true, plateKey);   % quiet

        if ~fitPrinted
            fprintf('\nRectangle clearance inside the plate ROI (px, must be >= 0):');
            for k = act, fprintf('  %s %.1f', rectNames{k}, info.rectFit(k)); end
            fprintf('\n\n');
            if any(info.rectFit(act) < 0)
                warning('An active rectangle sticks out of the plate ROI - shrink it or reduce g.inset in depth_geometry.m.');
            end
            fitPrinted = true;
        end

        ref = info.refMedian;  sg = info.refSigma;
        cutoff = ref - max(coldDelta, coldK * sg);   % background level
        if i == 1                                    % skip mask for cold spots that reach the plate edge
            skipMask = false(size(d));
            gr = stripSkipGrow;
            boxes = [geom.rectBox(setdiff(1:nR, act), :); geom.refExclude];
            for q = 1:size(boxes, 1)
                bx = boxes(q, :);
                skipMask(max(1, bx(3)-gr):min(end, bx(4)+gr), max(1, bx(1)-gr):min(end, bx(2)+gr)) = true;
            end
        end
        frac = zeros(1, 4);  drop = zeros(1, 4);
        for k = 1:4
            sm = info.mask & ~skipMask & (info.dist(:,:,k) < info.inset(k) + stripPx);   % strip just inside edge k
            s  = d(sm);
            frac(k) = mean(s < cutoff);
            drop(k) = ref - median(s, 'omitnan');
        end
        [worst, wk] = max(frac);
        [worstDrop, dk] = max(drop);

        % band from the reference median
        bandStr = 'OUT';
        for j = 1:size(bandRanges, 1)
            lo = bandRanges(j,1);  hi = bandRanges(j,2);
            inBand = ref >= lo && (ref < hi || (j == size(bandRanges,1) && ref <= hi));
            if inBand, bandStr = char(bandLabels(j)); break; end
        end

        bad = {};
        if worst > stripLimit, bad{end+1} = sprintf('%s edge', stripNames{wk}); end %#ok<SAGROW>
        dropLimit = max(medDrop, medDropK * sg);
        if worstDrop > dropLimit && ~(worst > stripLimit && wk == dk)
            bad{end+1} = sprintf('%s edge cooler', stripNames{dk}); %#ok<SAGROW>
        end
        for k = act                                  % active rectangles only
            if info.rect(k).pctValid < rectMinValid
                bad{end+1} = sprintf('%s mostly masked', rectNames{k}); %#ok<SAGROW>
            end
        end
        % frames in no band (hot start of a cooldown) skip the edge test
        if strcmp(bandStr, 'OUT')
            bad = {};  status = 'OUT (edge check skipped)';
        elseif isempty(bad)
            status = 'OK';
        else
            status = ['CHECK: ' strjoin(bad, ', ')];  nFlag = nFlag + 1;
        end
        bandOf(i) = bandStr;  usable(i) = isempty(bad) && ~strcmp(bandStr, 'OUT');
        rec.frame(i) = str2double(regexp(files(i).name, '\d+', 'match', 'once'));
        rec.ref(i) = ref;  rec.sigma(i) = sg;
        rec.edgeBgPct(i) = 100*worst;  rec.edgeDropC(i) = worstDrop;  rec.status(i) = string(status);
        for k = 1:nR
            rec.rMedian(i,k) = info.rect(k).median;  rec.rMean(i,k)  = info.rect(k).mean;
            rec.rStd(i,k)    = info.rect(k).std;     rec.rDT(i,k)    = info.rect(k).dT;
            rec.rValid(i,k)  = info.rect(k).pctValid;
        end

        fprintf('%-16s %-7.2f %-8s', files(i).name, ref, bandStr);
        for k = act, fprintf(' %-7.2f', rec.rDT(i,k)); end
        fprintf(' %-26s %s\n', sprintf('%s %.1f%%', stripNames{wk}, 100*worst), status);

        if showFlagged && ~isempty(bad)
            extract_depth_roi(d, true, files(i).name, true, plateKey);
        end
    end

    %3.Summary
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

    % dT per band
    fprintf('\n[%s] dT (degC), mean +/- sd over usable frames:\n', plateName);
    fprintf('  %-8s %3s', 'Band', 'n');
    for k = act, fprintf('  %-17s', rectNames{k}); end
    fprintf('\n');
    for lbl = bandLabels'
        sel = usable & (bandOf == lbl);
        if ~any(sel), continue; end
        fprintf('  %-8s %3d', lbl, nnz(sel));
        for k = act
            x = rec.rDT(sel, k);
            fprintf('  %7.2f +/- %-5.2f', mean(x), std(x));
        end
        fprintf('\n');
    end

    %4.Save
    T = table(rec.frame, bandOf, rec.ref, rec.sigma, rec.edgeBgPct, rec.edgeDropC, rec.status, ...
        'VariableNames', {'frame','band','refMedian','refSigma','edgeBgPct','edgeDropC','status'});
    for k = act
        nm = rectNames{k};
        T.([nm '_median']) = rec.rMedian(:,k);
        T.([nm '_mean'])   = rec.rMean(:,k);
        T.([nm '_std'])    = rec.rStd(:,k);
        T.([nm '_dT'])     = rec.rDT(:,k);
        T.([nm '_valid'])  = rec.rValid(:,k);
    end
    T = sortrows(T, {'band','frame'});
    csvPath = fullfile(statsFolder, ['depth_frame_bands_' plateName '.csv']);
    writetable(T, csvPath);
    fprintf('\nSaved band table: %s\n', csvPath);

    fprintf('%d of %d frames flagged in %s.\n', nFlag, numel(files), plateName);

    %5.Band figures (the edge test only catches a plate shift of about 10 px or more)
    if showPerBand
        for lbl = bandLabels'
            idx = find(bandOf == lbl, 1);
            if isempty(idx), continue; end
            dd = load_thermal_matrix(fullfile(folder, files(idx).name));
            extract_depth_roi(dd, true, sprintf('%s  %s  band %s', plateName, files(idx).name, lbl), true, plateKey);
        end
    end
end