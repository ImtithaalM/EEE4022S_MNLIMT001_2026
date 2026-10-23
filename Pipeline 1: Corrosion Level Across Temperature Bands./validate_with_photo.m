% VALIDATE_WITH_PHOTO  Validation test: do the thermal anomalies sit on corrosion that is VISIBLE in the plate photo?
%
% Plates: light (localised patch), FULL (uniform corrosion) and HALF (partly corroded plate).  The photo is only the answer key.
% thermal_anomaly_map already decides what is anomalous; here that map is laid over the photo and checked:
%   1) Registration: 4 clicks on the photo (the plate corners, in the orientation of the THERMAL image) line the
%      photo up with the thermal ROI exactly (the plate outline and ROI insets from extract_plate_roi are used,
%      so no tape or extra landmarks are needed). Clicks are saved per plate and reused.
%   2) Rust-coloured photo pixels are found. If stats\rust_protos_<plate>.mat exists (from tune_rust) your tuned
%      nearest-example classifier is used (rust_colour_rel); otherwise a simple Lab colour rule (brown/orange =
%      positive a* and b*), which works for vivid rust such as the light plate.  STRICT and LOOSE levels as before.
%   3) Over every frame of the plate in anomaly_masks%maskTag%.mat:
%        - pixel scores: mean local rust density at hot / cold / anomalous / non-anomalous pixels, LIFT (> 1 = the
%          anomalies sit on more rust than the plate average) and Spearman rho (rust density vs thermal residual);
%        - overlap of the LARGEST joined hot region (the localised-corrosion indicator) with rust:
%          precision = share of that region lying on rust, recall = share of rust area covered, Dice.
%      A tolerance (tolPx) widens the rust area a little so a few pixels of registration error do not count against it.
%   Three rust levels: STRICT (clear rust), LOOSE (also duller rust), FAINT (lighter corrosion around a main patch).
%   Expect: LIGHT -> large hot region on the rust patch (high precision / recall / lift).  FULL -> almost no
%   anomalies (uniform), so overlap is not defined; the check there is photo rust % vs thermal (rho, rust density).
%
% BAND GROUPS (22-24, 24-26, 28-32, 38-42 and 48-52): the photo and the 4 corner clicks are per PLATE, so each plate is registered once and then
% scored separately for every band group found in the masks (one results row and one figure per plate and group).
% The 48-52 group is the POOLED band made by thermal_anomaly_map (48-50 and 50-52 merged, HALF 1573 excluded there).
% 22-24 and 24-26 are the two room-temperature groups of thermal_anomaly_map. FULL has no frames in 22-24, light only 2, and HALF only 2 in 24-26,
% so those rows rest on very few frames. The corner clicks are always taken on the 28-32 reference frame when it exists (saved clicks are reused).
% Run thermal_anomaly_map first (saves stats\anomaly_masks_withTape.mat). Needs the Image Processing Toolbox.
% Put 10mm_light.jpg, 10mm_full.jpg and 10mm_half.jpg (HALF) in the photo folder (DIGITAL).
% HALF: first run asks for 4 corner clicks (saved to photo_registration_points_HALF.mat). Run thermal_anomaly_map with HALF first.
% Output: figures and a results table (one row per plate and band group) in results_photo_check\

% ---- settings ----
%root = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\matlab pipeline';
root        = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\FINAL PIPELINES\corr pp';
statsFolder = fullfile(root, 'stats');
outFolder   = fullfile(root, 'results_photo_check');
if ~isfolder(outFolder), mkdir(outFolder); end
photoDir    = fullfile(root, 'DIGITAL');
plates = struct('name',     {'light',           'FULL',           'HALF'}, ...
                'photo',    {'10mm_light.jpg',  '10mm.png',  '10mm_half.jpg'}, ...
                'refFrame', {{1898 1883 1853 1810 1786}, {0 1770 1756 1749 1737}, {1719 1715 1676 1620 1580}}, ...   % reference frame per band group (order of groupTags; 0 = plate has no frames in that group)
                'bias',     {[0 -1],            [-2 -4],          [0 -1]}, ...     % [strict loose] bias for the tuned classifier (higher = stricter)
                'method',   {'lab',             'lab',            'dark'}, ...
                'map',      {'local',           'local',          'ref'}, ...      % thermal map used: 'local' = local-residual anomalies (hot + cold), 'ref' = steel-referenced hot area (for a plate with a large corroded share)
                'warpBlur', {1.5,               1.5,              0.6}, ...        % blur (photo px) applied to the rust maps before warping; smaller keeps finer rust
                'warpThr',  {0.5,               0.5,              0.3});           % a warped pixel counts as rust above this (lower keeps finer rust)        % rust rule when there is no tuned classifier: 'lab' = a*/b* colour, 'dark' = darker than the plate's bare steel
darkThr     = [-6 0 4];        % 'dark' rule: strict / loose / faint = Otsu lightness threshold + these (L units); rust must also have a* above [2 1 0]
darkA       = [2 1 0];
labThr      = [8 12 6 9 3 6];  % Lab fallback: strict a*, b*, loose a*, b*, FAINT a*, b* (rust = a* and b* above these)
maskTag     = '_withTape';     % which anomaly_masks file
groupTags   = {'22-24', '24-26', '28-32', '38-42', '48-52'};   % band groups, in the order of the refFrame lists above
useGroups   = groupTags;        % validate ALL band groups (one figure + one results row per plate and band group); e.g. {'38-42'} for one band only
regGroup    = '38-42';         % band group whose reference frame is shown for the corner clicks (must be in useGroups)
photoMaxSide = 1000;           % photo is shrunk to this on its long side for the colour work
densSigma   = 8;               % rust becomes a local density (Gaussian, ROI pixels) so registration error is forgiven
tolPx       = 6;               % binary overlap tolerance (ROI pixels)
closeRad    = 4;  bigMin = 500;   % same joining rule as thermal_anomaly_map (largest region)
redoPoints  = false;           % true = click the photo corners again (set back to false after you have clicked)
nLarge      = 5;               % frames used for the large-scale check (raw temperature vs rust, needs the xlsx re-read)
sigmaLarge  = 20;              % px: smoothing scale of that check
% plate outline used by extract_plate_roi (MUST match its OUT_* values)
OUT_TOP = 59;  OUT_BOTTOM = 434;  OUT_LEFT = 56;  OUT_RIGHT = 579;

% ---- thermal anomaly data ----
A = load(fullfile(statsFolder, ['anomaly_masks' maskTag '.mat']));
jobs = A.jobs;
if size(jobs, 2) < 5, jobs(:,5) = {'28-32'}; end      % masks saved before the band groups existed
% The saved masks hold the frame folders of the computer they were made on (thermal_anomaly_map). Rebuild them from 'root'
% so this script runs on any computer.
for j = 1:size(jobs, 1)
    nm = jobs{j,1};  if strcmp(nm, 'bare (control)'), nm = 'bare'; end
    jobs{j,2} = fullfile(root, nm);
end
summ = {};

for pi_ = 1:numel(plates)
    pl = plates(pi_);
    jP = find(strcmp(jobs(:,1), pl.name));
    if isempty(jP), warning('No %s frames in anomaly_masks%s.mat - skipped.', pl.name, maskTag); continue; end
    gAvail = groupTags(ismember(groupTags, jobs(jP,5)) & ismember(groupTags, useGroups));          % band groups that have frames of this plate
    gReg = regGroup;  if ~ismember(gReg, gAvail), gReg = gAvail{1}; end
    jG1  = jP(strcmp(jobs(jP,5), gReg));                           % registration: reference frame of the registration group (one photo per plate)
    jRef = jG1(find(cell2mat(jobs(jG1,3)) == pl.refFrame{find(strcmp(groupTags, gReg), 1)}, 1));
    if isempty(jRef), jRef = jG1(1); end
    refFrame = jobs{jRef,3};
    [H, W] = size(A.resS{jRef});

    % ROI position inside the thermal frame
    d0 = load_thermal_matrix(fullfile(jobs{jRef,2}, sprintf('IV_%05d.xlsx', refFrame)));
    [~, info] = extract_plate_roi(d0, false, '', true);
    rows = info.rows;  cols = info.cols;
    if numel(rows) ~= H || numel(cols) ~= W, error('ROI size does not match the saved masks: rerun thermal_anomaly_map.'); end

    % ---- photo and corner clicks ----
    img = imread(fullfile(photoDir, pl.photo));
    ptsFile = fullfile(statsFolder, ['photo_registration_points_' pl.name '.mat']);
    doClick = true;
    if isfile(ptsFile) && ~redoPoints
        P = load(ptsFile);
        if isfield(P, 'photoPts') && isequal(size(P.photoPts), [4 2]), photoPts = P.photoPts;  doClick = false; end
    end
    if doClick       % no usable saved clicks for this plate (or redoPoints = true)
        labels = {'TOP-LEFT', 'TOP-RIGHT', 'BOTTOM-RIGHT', 'BOTTOM-LEFT'};
        fig = figure('Name', ['Click plate corners: ' pl.name], 'Position', [50 50 1400 800]);
        subplot(1,2,2);  h = imagesc(double(A.resS{jRef}), [-5 5]);  set(h, 'AlphaData', ~isnan(double(A.resS{jRef})));
        set(gca, 'Color', [0.85 0.85 0.85]);  axis image;  title(sprintf('Thermal %s IV\\_%05d: top-left of THIS image is TOP-LEFT', pl.name, refFrame));
        axP = subplot(1,2,1);  imshow(img, 'Parent', axP);  hold(axP, 'on');
        photoPts = zeros(4, 2);
        for k = 1:4
            okClick = false;
            while ~okClick
                title(axP, sprintf('Click %d of 4 ON THIS PHOTO (left): the plate corner that is %s in the THERMAL image (right)', k, labels{k}));
                ginput(1);
                cp = axP.CurrentPoint(1, 1:2);  x = cp(1);  y = cp(2);
                xl = xlim(axP);  yl = ylim(axP);
                okClick = x >= xl(1) && x <= xl(2) && y >= yl(1) && y <= yl(2);
                if ~okClick, title(axP, 'That click was NOT on the photo - click the PHOTO on the left'); pause(1.2); end
            end
            photoPts(k, :) = [x y];
            plot(axP, x, y, 'r+', 'MarkerSize', 14, 'LineWidth', 2);  text(axP, x + 10, y, num2str(k), 'Color', 'r', 'FontSize', 14);
        end
        close(fig);
        save(ptsFile, 'photoPts');
    end

    % ---- shrink photo, rust colour ----
    s = photoMaxSide / max(size(img, 1), size(img, 2));
    imgS = imresize(img, s);
    ptsS = photoPts * s;
    [hS, wS, ~] = size(imgS);
    plateS = imerode(poly2mask(ptsS(:,1), ptsS(:,2), hS, wS), strel('disk', 6));
    protoFile = fullfile(statsFolder, ['rust_protos_' pl.name '.mat']);
    if isfile(protoFile)
        PP = load(protoFile);
        [rsS, rlS] = rust_colour_rel(imgS, plateS, PP.steelP, PP.rustP, pl.bias(1), pl.bias(2));
        [~, rfS]   = rust_colour_rel(imgS, plateS, PP.steelP, PP.rustP, pl.bias(2), pl.bias(2) - 1);   % FAINT = one step looser than loose
        method = 'tuned classifier (tune_rust)';
    else
        if strcmp(pl.method, 'dark')
            [rsS, rlS, rfS] = rust_colour_dark(imgS, plateS, darkThr, darkA);
            method = 'lightness rule: rust = darker than the bare steel (no rust_protos file)';
        else
            [rsS, rlS, rfS] = rust_colour_lab(imgS, plateS, labThr);
            method = 'Lab colour rule (no rust_protos file)';
        end
    end
    photoRustStrict = 100 * nnz(rsS) / nnz(plateS);
    photoRustLoose  = 100 * nnz(rlS) / nnz(plateS);
    photoRustFaint  = 100 * nnz(rfS) / nnz(plateS);

    % ---- register to the thermal ROI (outline corners in thermal-frame pixels) ----
    outline = [OUT_LEFT OUT_TOP; OUT_RIGHT OUT_TOP; OUT_RIGHT OUT_BOTTOM; OUT_LEFT OUT_BOTTOM];
    tform  = fitgeotrans(ptsS, outline, 'projective');
    outRef = imref2d([H W], [cols(1)-0.5 cols(end)+0.5], [rows(1)-0.5 rows(end)+0.5]);
    photoW     = imwarp(imgS, tform, 'OutputView', outRef);
    photoValid = imwarp(plateS, tform, 'OutputView', outRef) > 0.5;
    rustStrict = imwarp(imgaussfilt(double(rsS), pl.warpBlur), tform, 'OutputView', outRef) > pl.warpThr;
    rustLoose  = imwarp(imgaussfilt(double(rlS), pl.warpBlur), tform, 'OutputView', outRef) > pl.warpThr;
    rustFaint  = imwarp(imgaussfilt(double(rfS), pl.warpBlur), tform, 'OutputView', outRef) > pl.warpThr;
    rustMaps = {rustStrict, rustLoose, rustFaint};  mapNames = {'STRICT', 'LOOSE', 'FAINT'};
    vD = double(photoValid);
    tolMaps = cellfun(@(m) imgaussfilt(double(m) .* vD, densSigma) ./ max(imgaussfilt(vD, densSigma), eps), rustMaps, 'UniformOutput', false);

    for gq = 1:numel(gAvail)
        gt = gAvail{gq};  gi = find(strcmp(groupTags, gt), 1);
        jG = jP(strcmp(jobs(jP,5), gt));                          % this plate's frames in this band group
        jRefG = jG(find(cell2mat(jobs(jG,3)) == pl.refFrame{gi}, 1));
        if isempty(jRefG), jRefG = jG(1); end
        refFrameG = jobs{jRefG,3};
        % ---- scores over every frame of this plate ----
        nF = numel(jG);
        sc  = nan(nF, 3, 7);     % frame x map x [base pHot pCold pAny pNon lift rho]
        ov  = nan(nF, 3, 3);     % frame x rust level x [precision recall dice] of the largest hot region
        ovA = nan(nF, 1);        % size of the largest hot region, % of the plate
        thermStat = nan(nF, 3);  % [hot% cold% largest-region%] of the ROI, thermal side of the photo comparison
        rhoL = nan(nF, 1);       % large-scale Spearman rho: smoothed RAW temperature vs smoothed rust density
        idxLarge = unique(round(linspace(1, nF, min(nLarge, nF))));
        for q = 1:nF
            j = jG(q);
            r = double(A.resS{j});
            valid = ~isnan(r) & photoValid;
            if nnz(valid) < 50, fprintf('Frame %d skipped: only %d valid px\n', jobs{j,3}, nnz(valid)); continue, end
            if strcmp(pl.map, 'ref')
                hot = A.hotRefM{j} & valid;  cold = false(size(hot));      % steel-referenced map: hot area only
            else
                hot = A.hotM{j} & valid;  cold = A.coldM{j} & valid;
            end
            anyA = hot | cold;
            thermStat(q, :) = [100*nnz(hot)/nnz(valid), 100*nnz(cold)/nnz(valid), 0];
            for m = 1:3
                rt = tolMaps{m};
                base = mean(rt(valid));
                pAny = mean(rt(anyA));
                sc(q, m, :) = [base, mean(rt(hot)), mean(rt(cold)), pAny, mean(rt(valid & ~anyA)), pAny / base, ...
                               corr(rt(valid), r(valid), 'Type', 'Spearman')];
            end
            bigM = big_regions(hot, closeRad, bigMin);
            if any(bigM(:))
                for m = 1:3
                    rb = rustMaps{m} & photoValid;
                    prec = nnz(bigM & imdilate(rb, strel('disk', tolPx))) / nnz(bigM);
                    rec  = nnz(rb & imdilate(bigM, strel('disk', tolPx))) / max(nnz(rb), 1);
                    ov(q, m, :) = [prec, rec, 2*prec*rec / max(prec + rec, eps)];
                end
                ovA(q) = 100 * nnz(bigM) / nnz(valid);
                thermStat(q, 3) = ovA(q);
            end
            if ismember(q, idxLarge)      % large-scale check on the raw (not high-passed) temperature
                [rawROI, ~] = extract_plate_roi(load_thermal_matrix(fullfile(jobs{j,2}, sprintf('IV_%05d.xlsx', jobs{j,3}))), false, '', true);
                vv = ~isnan(rawROI) & photoValid;
                tS = nan_smooth(rawROI, vv, sigmaLarge);
                rS = nan_smooth(double(rustLoose), photoValid, sigmaLarge);
                rhoL(q) = corr(tS(vv), rS(vv), 'Type', 'Spearman');
            end
        end
        M  = squeeze(mean(sc, 1, 'omitnan'));      % 3 x 7
        Ov = squeeze(mean(ov, 1, 'omitnan'));      % 3 x 3  (rust level x [precision recall dice])
        nBig = nnz(~isnan(ov(:,1,1)));
        bigPct = mean(ovA, 'omitnan');
        Th = mean(thermStat, 1, 'omitnan');                 % thermal: hot%, cold%, largest region %
        rhoLarge = mean(rhoL, 'omitnan');
        rustPct = [100*nnz(rustStrict & photoValid), 100*nnz(rustLoose & photoValid), 100*nnz(rustFaint & photoValid)] / nnz(photoValid);   % photo rust, % of the ROI

        fprintf('\n=== %s, band group %s: do the thermal anomalies sit on rust visible in the photo? (%d frames) ===\n', upper(pl.name), gt, nF);
        fprintf('Rust colour method: %s.  Photo rust area: strict %.1f %%, loose %.1f %%, faint %.1f %% of the plate.\n', method, photoRustStrict, photoRustLoose, photoRustFaint);
        fprintf('%-8s %-8s %-12s %-13s %-13s %-13s %-6s %-6s\n', 'Rule', 'Rust %', 'P(rust|hot)', 'P(rust|cold)', 'P(rust|any)', 'P(rust|none)', 'LIFT', 'Rho');
        for m = 1:3
            fprintf('%-8s %-8.1f %-12.1f %-13.1f <strong>%-13.1f</strong> %-13.1f <strong>%-6.2f</strong> %-6.2f\n', mapNames{m}, ...
                100*M(m,1), 100*M(m,2), 100*M(m,3), 100*M(m,4), 100*M(m,5), M(m,6), M(m,7));
        end
        if nBig > 0
            fprintf('Largest hot region (%d of %d frames have one, %.1f %% of the plate) vs each rust level:\n', nBig, nF, bigPct);
            fprintf('  %-8s precision %.2f   recall %.2f   Dice %.2f\n', 'STRICT', Ov(1,1), Ov(1,2), Ov(1,3));
            fprintf('  %-8s <strong>precision %.2f</strong>   recall %.2f   Dice %.2f\n', 'LOOSE', Ov(2,1), Ov(2,2), Ov(2,3));
            fprintf('  %-8s precision %.2f   recall %.2f   Dice %.2f   (faint rust: lighter corrosion around the main patch)\n', 'FAINT', Ov(3,1), Ov(3,2), Ov(3,3));
        else
            fprintf('No large hot region in any frame: overlap not defined (uniform plate, as expected for FULL).\n');
        end
        fprintf('\nPHOTO RUST vs THERMAL ANOMALY (%% of the ROI, mean over %d frames):\n', nF);
        fprintf('  Photo rust    strict %5.1f   loose %5.1f   faint %5.1f\n', rustPct(1), rustPct(2), rustPct(3));
        fprintf('  Thermal       hot %5.1f   cold %5.1f   anomalous total %5.1f   largest hot region <strong>%5.1f</strong>\n', Th(1), Th(2), Th(1) + Th(2), Th(3));
        fprintf('  Thermal anomaly / photo rust (faint) = %.2f,  (loose) = %.2f\n', (Th(1) + Th(2)) / max(rustPct(3), eps), (Th(1) + Th(2)) / max(rustPct(2), eps));
        fprintf('  Large-scale check (raw temperature vs rust, %d px smoothing, %d frames): Spearman rho = %.2f\n', sigmaLarge, numel(idxLarge), rhoLarge);
        fprintf('Densities are mean local rust density (%%) at those pixels. LIFT > 1 = anomalies sit on more rust than average.\n');
        fprintf('Rho = Spearman correlation of rust density with the thermal residual (+ = rusty areas are warmer).\n');
        summ(end+1, :) = {pl.name, gt, nF, rustPct(1), rustPct(2), rustPct(3), Th(1), Th(2), Th(3), rhoLarge, 100*M(2,1), 100*M(2,4), M(2,6), M(2,7), nBig, bigPct, ...
            Ov(2,1), Ov(2,2), Ov(2,3), Ov(3,1), Ov(3,2), Ov(3,3)}; %#ok<SAGROW>

        % ---- figure for the reference frame ----
        r = double(A.resS{jRefG});
        valid = ~isnan(r) & photoValid;
        if strcmp(pl.map, 'ref')
            hot = A.hotRefM{jRefG} & valid;  cold = false(size(hot));
        else
            hot = A.hotM{jRefG} & valid;  cold = A.coldM{jRefG} & valid;
        end
        bigM = big_regions(hot, closeRad, bigMin);
        f = figure('Name', ['Validation ' pl.name ' ' gt], 'Position', [40 40 1500 460]);
        tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        ax1 = nexttile(tl);  imshow(photoW, 'Parent', ax1);
        title(ax1, {sprintf('%s: photo', pl.name), 'registered to the ROI'}, 'FontSize', 10);
        ax2 = nexttile(tl);  imshow(labeloverlay(photoW, uint8(rustLoose & photoValid) + uint8(rustFaint & photoValid), 'Colormap', [1 0.7 1; 1 0 1], 'Transparency', 0.45), 'Parent', ax2);
        title(ax2, {'Photo rust (% of ROI):', sprintf('loose %.1f (magenta), faint %.1f (pale)', rustPct(2), rustPct(3))}, 'FontSize', 10);
        ax3 = nexttile(tl);  h = imagesc(ax3, r, [-5 5]);  set(h, 'AlphaData', ~isnan(r));  set(ax3, 'Color', [0.85 0.85 0.85]);  axis(ax3, 'image');
        cb = colorbar(ax3, 'southoutside');  cb.Label.String = 'residual (\circC)';
        title(ax3, {sprintf('Thermal residual, IV\\_%05d:', refFrameG), sprintf('anomalous %.1f %% of ROI, largest region %.1f %%', Th(1) + Th(2), Th(3))}, 'FontSize', 10);
        for ax = findall(f, 'Type', 'axes')', try, ax.Toolbar = []; catch, end; end
        stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        try
            exportgraphics(f, fullfile(outFolder, ['photo_validation_' pl.name '_' gt maskTag '_' stamp '.png']), 'Resolution', 300);
        catch ME
            fprintf('Figure for %s NOT saved (%s).\n', pl.name, ME.message);
        end
    end
end

% ---- results table ----
if ~isempty(summ)
    T = cell2table(summ, 'VariableNames', {'Plate','Group','Frames','RoiRustStrictPct','RoiRustLoosePct','RoiRustFaintPct','ThermHotPct','ThermColdPct','ThermLargestRegionPct','RhoLargeScale','RustDensityPct','PRustGivenAnomalyPct', ...
        'Lift','Rho','FramesWithBigRegion','BigRegionPctOfPlate','PrecisionLoose','RecallLoose','DiceLoose','PrecisionFaint','RecallFaint','DiceFaint'});
    writetable(T, fullfile(outFolder, ['photo_validation_results' maskTag '.csv']));
    fprintf('\nSaved: %s\n', fullfile(outFolder, ['photo_validation_results' maskTag '.csv']));
end

% ---------------------------------------------------------------------
function [rs, rl, rf] = rust_colour_lab(rgb, useMask, thr)
% Rust = brown / orange: Lab a* and b* both positive and above the thresholds (strict / loose / faint).
lab = rgb2lab(im2double(rgb));
Af = imgaussfilt(lab(:,:,2), 1.2);
Bf = imgaussfilt(lab(:,:,3), 1.2);
rs = (Af > thr(1)) & (Bf > thr(2)) & useMask;
rl = (Af > thr(3)) & (Bf > thr(4)) & useMask;
rf = (Af > thr(5)) & (Bf > thr(6)) & useMask;
end

function [rs, rl, rf] = rust_colour_dark(rgb, useMask, offs, aMin)
% Rust = darker than the plate's bare steel. The lightness threshold is the Otsu split of L* over the plate, so it adapts to the
% photo. Needed when the rust is dark brown with little colour (a*, b* then look like the steel). Rust must also lean red (a* > aMin).
% FINE rust (small pits and specks, thin streaks) is added at the loose / faint levels: a pixel clearly darker than the steel
% around it (local median over about 80 px) counts too, even if it is not below the global threshold.
lab = rgb2lab(im2double(rgb));
Lf = imgaussfilt(lab(:,:,1), 0.8);
Af = imgaussfilt(lab(:,:,2), 0.8);
t = 100 * graythresh(Lf(useMask) / 100);
Ls = imresize(Lf, 0.25);
Lbg = imresize(medfilt2(Ls, [21 21], 'symmetric'), size(Lf), 'bilinear');      % local steel level
fineL = (Lf < Lbg - 6) & (Af > -1) & (Lf < t + 14);        % clearly darker than its surroundings (loose)
fineF = (Lf < Lbg - 4) & (Af > -1) & (Lf < t + 14);        % a little less strict (faint)
rs = (Lf < t + offs(1)) & (Af > aMin(1)) & useMask;
rl = (((Lf < t + offs(2)) & (Af > aMin(2))) | fineL) & useMask;
rf = (((Lf < t + offs(3)) & (Af > aMin(3))) | fineF | ((Lf < t + 10) & (Af > -1))) & useMask;
end

function bigM = big_regions(mask, rad, minA)
% Anomalous pixels belonging to joined regions of at least minA pixels (same joining as thermal_anomaly_map).
cc = bwconncomp(imclose(mask, strel('disk', rad)));
a = cellfun(@numel, cc.PixelIdxList);
bigM = false(size(mask));
for i = find(a >= minA)
    bigM(cc.PixelIdxList{i}) = true;
end
bigM = bigM & mask;
end

function y = nan_smooth(x, valid, sg)
% Gaussian smoothing that ignores pixels outside 'valid' (NaN-aware).
z = x;  z(~valid) = 0;
w = imgaussfilt(double(valid), sg);
y = imgaussfilt(z, sg) ./ max(w, eps);
y(~valid) = NaN;
end