% VALIDATE_WITH_PHOTO  Checks whether the thermal anomalies from thermal_anomaly_map sit on corrosion visible in the plate photo (light, FULL, HALF).
% Registers each photo to the thermal ROI with 4 corner clicks (saved per plate), finds rust colour, then scores every frame per band group (rust lift, Spearman rho, precision / recall / Dice of the largest hot region).
% Run thermal_anomaly_map first; photos go in DIGITAL. Output: figures and photo_validation_results<tag>.csv in results_photo_check\

% ---- settings ----
%root = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\matlab pipeline';
root        = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\matlab pipeline';
statsFolder = fullfile(root, 'stats');
outFolder   = fullfile(root, 'results_photo_check');
if ~isfolder(outFolder), mkdir(outFolder); end
photoDir    = fullfile(root, 'DIGITAL');
% refFrame = reference frame per band group (order of groupTags; 0 = no frames); bias = [strict loose] for the tuned classifier;
% method = rust rule without a tuned classifier ('lab' or 'dark'); map = 'local' (hot + cold) or 'ref' (steel-referenced hot area);
% warpBlur / warpThr = blur and threshold applied to the rust maps when warping (lower keeps finer rust)
plates = struct('name',     {'light',           'FULL',           'HALF'}, ...
                'photo',    {'10mm_light.jpg',  '10mm_full.jpg',  '10mm_half.jpg'}, ...
                'refFrame', {{1898 1883 1853 1810 1786}, {0 1770 1756 1749 1737}, {1719 1715 1680 1620 1580}}, ...
                'bias',     {[0 -1],            [-2 -4],          [0 -1]}, ...
                'method',   {'lab',             'lab',            'dark'}, ...
                'map',      {'local',           'local',          'ref'}, ...
                'warpBlur', {1.5,               1.5,              0.6}, ...
                'warpThr',  {0.5,               0.5,              0.3});
darkThr     = [-6 0 4];        % 'dark' rule: Otsu lightness threshold + these (strict / loose / faint)
darkA       = [2 1 0];         % 'dark' rule: minimum a* (strict / loose / faint)
labThr      = [8 12 6 9 3 6];  % 'lab' rule: a*, b* thresholds (strict, loose, faint)
maskTag     = '_withTape';     % which anomaly_masks file
groupTags   = {'22-24', '24-26', '28-32', '38-42', '48-52'};   % band groups, in the order of the refFrame lists
regGroup    = '28-32';         % group whose reference frame is shown for the corner clicks
photoMaxSide = 1000;           % photo long side (px)
densSigma   = 8;               % rust density smoothing (ROI px)
tolPx       = 6;               % binary overlap tolerance (ROI px)
closeRad    = 4;  bigMin = 500;   % joining rule for the largest region (as thermal_anomaly_map)
redoPoints  = false;           % true = click the photo corners again
nLarge      = 5;               % frames used for the large-scale check
sigmaLarge  = 20;              % px: smoothing scale of the large-scale check
% plate outline (must match extract_plate_roi)
OUT_TOP = 59;  OUT_BOTTOM = 434;  OUT_LEFT = 56;  OUT_RIGHT = 579;

%1.Thermal anomaly data
A = load(fullfile(statsFolder, ['anomaly_masks' maskTag '.mat']));
jobs = A.jobs;
if size(jobs, 2) < 5, jobs(:,5) = {'28-32'}; end
for j = 1:size(jobs, 1)
    nm = jobs{j,1};  if strcmp(nm, 'bare (control)'), nm = 'bare'; end
    jobs{j,2} = fullfile(root, nm);
end
summ = {};

for pi_ = 1:numel(plates)
    pl = plates(pi_);
    jP = find(strcmp(jobs(:,1), pl.name));
    if isempty(jP), warning('No %s frames in anomaly_masks%s.mat - skipped.', pl.name, maskTag); continue; end
    gAvail = groupTags(ismember(groupTags, jobs(jP,5)));
    gReg = regGroup;  if ~ismember(gReg, gAvail), gReg = gAvail{1}; end
    jG1  = jP(strcmp(jobs(jP,5), gReg));
    jRef = jG1(find(cell2mat(jobs(jG1,3)) == pl.refFrame{find(strcmp(groupTags, gReg), 1)}, 1));
    if isempty(jRef), jRef = jG1(1); end
    refFrame = jobs{jRef,3};
    [H, W] = size(A.resS{jRef});

    % ROI position inside the thermal frame
    d0 = load_thermal_matrix(fullfile(jobs{jRef,2}, sprintf('IV_%05d.xlsx', refFrame)));
    [~, info] = extract_plate_roi(d0, false, '', true);
    rows = info.rows;  cols = info.cols;
    if numel(rows) ~= H || numel(cols) ~= W, error('ROI size does not match the saved masks: rerun thermal_anomaly_map.'); end

    %2.Photo and corner clicks
    img = imread(fullfile(photoDir, pl.photo));
    ptsFile = fullfile(statsFolder, ['photo_registration_points_' pl.name '.mat']);
    doClick = true;
    if isfile(ptsFile) && ~redoPoints
        P = load(ptsFile);
        if isfield(P, 'photoPts') && isequal(size(P.photoPts), [4 2]), photoPts = P.photoPts;  doClick = false; end
    end
    if doClick
        labels = {'TOP-LEFT', 'TOP-RIGHT', 'BOTTOM-RIGHT', 'BOTTOM-LEFT'};
        fig = figure('Name', ['Click plate corners: ' pl.name], 'Position', [50 50 1400 800]);
        subplot(1,2,2);  h = imagesc(double(A.resS{jRef}), [-5 5]);  set(h, 'AlphaData', ~isnan(double(A.resS{jRef})));
        set(gca, 'Color', [0.85 0.85 0.85]);  axis image;  title(sprintf('Thermal %s IV\\_%05d: top-left of THIS image is TOP-LEFT', pl.name, refFrame));
        subplot(1,2,1);  imshow(img);  hold on;
        photoPts = zeros(4, 2);
        for k = 1:4
            title(sprintf('Click %d of 4: the plate corner that is %s in the THERMAL image (right)', k, labels{k}));
            [x, y] = ginput(1);
            photoPts(k, :) = [x y];
            plot(x, y, 'r+', 'MarkerSize', 14, 'LineWidth', 2);  text(x + 10, y, num2str(k), 'Color', 'r', 'FontSize', 14);
        end
        close(fig);
        save(ptsFile, 'photoPts');
    end

    %3.Rust colour
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

    %4.Register to the thermal ROI
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
        jG = jP(strcmp(jobs(jP,5), gt));
        jRefG = jG(find(cell2mat(jobs(jG,3)) == pl.refFrame{gi}, 1));
        if isempty(jRefG), jRefG = jG(1); end
        refFrameG = jobs{jRefG,3};

        %5.Scores over every frame
        nF = numel(jG);
        sc  = nan(nF, 3, 7);     % frame x map x [base pHot pCold pAny pNon lift rho]
        ov  = nan(nF, 3, 3);     % frame x rust level x [precision recall dice] of the largest hot region
        ovA = nan(nF, 1);        % largest hot region, % of the plate
        thermStat = nan(nF, 3);  % [hot% cold% largest-region%]
        rhoL = nan(nF, 1);       % large-scale rho: smoothed raw temperature vs smoothed rust density
        idxLarge = unique(round(linspace(1, nF, min(nLarge, nF))));
        for q = 1:nF
            j = jG(q);
            r = double(A.resS{j});
            valid = ~isnan(r) & photoValid;
            if strcmp(pl.map, 'ref')
                hot = A.hotRefM{j} & valid;  cold = false(size(hot));
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
            if ismember(q, idxLarge)
                [rawROI, ~] = extract_plate_roi(load_thermal_matrix(fullfile(jobs{j,2}, sprintf('IV_%05d.xlsx', jobs{j,3}))), false, '', true);
                vv = ~isnan(rawROI) & photoValid;
                tS = nan_smooth(rawROI, vv, sigmaLarge);
                rS = nan_smooth(double(rustLoose), photoValid, sigmaLarge);
                rhoL(q) = corr(tS(vv), rS(vv), 'Type', 'Spearman');
            end
        end
        M  = squeeze(mean(sc, 1, 'omitnan'));      % 3 x 7
        Ov = squeeze(mean(ov, 1, 'omitnan'));      % 3 x 3
        nBig = nnz(~isnan(ov(:,1,1)));
        bigPct = mean(ovA, 'omitnan');
        Th = mean(thermStat, 1, 'omitnan');
        rhoLarge = mean(rhoL, 'omitnan');
        rustPct = [100*nnz(rustStrict & photoValid), 100*nnz(rustLoose & photoValid), 100*nnz(rustFaint & photoValid)] / nnz(photoValid);

        %6.Results
        fprintf('\n=== %s, band group %s (%d frames) ===\n', upper(pl.name), gt, nF);
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
            fprintf('  %-8s precision %.2f   recall %.2f   Dice %.2f\n', 'FAINT', Ov(3,1), Ov(3,2), Ov(3,3));
        else
            fprintf('No large hot region in any frame: overlap not defined.\n');
        end
        fprintf('\nPHOTO RUST vs THERMAL ANOMALY (%% of the ROI, mean over %d frames):\n', nF);
        fprintf('  Photo rust    strict %5.1f   loose %5.1f   faint %5.1f\n', rustPct(1), rustPct(2), rustPct(3));
        fprintf('  Thermal       hot %5.1f   cold %5.1f   anomalous total %5.1f   largest hot region <strong>%5.1f</strong>\n', Th(1), Th(2), Th(1) + Th(2), Th(3));
        fprintf('  Thermal anomaly / photo rust (faint) = %.2f,  (loose) = %.2f\n', (Th(1) + Th(2)) / max(rustPct(3), eps), (Th(1) + Th(2)) / max(rustPct(2), eps));
        fprintf('  Large-scale check (raw temperature vs rust, %d px smoothing, %d frames): Spearman rho = %.2f\n', sigmaLarge, numel(idxLarge), rhoLarge);
        summ(end+1, :) = {pl.name, gt, nF, rustPct(1), rustPct(2), rustPct(3), Th(1), Th(2), Th(3), rhoLarge, 100*M(2,1), 100*M(2,4), M(2,6), M(2,7), nBig, bigPct, ...
            Ov(2,1), Ov(2,2), Ov(2,3), Ov(3,1), Ov(3,2), Ov(3,3)}; %#ok<SAGROW>

        %7.Figure for the reference frame
        r = double(A.resS{jRefG});
        valid = ~isnan(r) & photoValid;
        if strcmp(pl.map, 'ref')
            hot = A.hotRefM{jRefG} & valid;  cold = false(size(hot));
        else
            hot = A.hotM{jRefG} & valid;  cold = A.coldM{jRefG} & valid;
        end
        bigM = big_regions(hot, closeRad, bigMin);
        f = figure('Name', ['Validation ' pl.name ' ' gt], 'Position', [40 40 1400 800]);
        subplot(2,2,1);  imshow(photoW);  title(sprintf('%s photo warped to the thermal ROI', pl.name));
        subplot(2,2,2);  imshow(labeloverlay(photoW, uint8(rustLoose & photoValid) + uint8(rustFaint & photoValid), 'Colormap', [1 0.7 1; 1 0 1], 'Transparency', 0.45));
        title(sprintf('Photo rust (%% of ROI): loose %.1f (magenta), faint %.1f (pale)', rustPct(2), rustPct(3)));
        subplot(2,2,3);  h = imagesc(r, [-5 5]);  set(h, 'AlphaData', ~isnan(r));  set(gca, 'Color', [0.85 0.85 0.85]);  axis image;  colorbar;
        title(sprintf('Thermal residual, IV\\_%05d: anomalous %.1f %% of ROI, largest region %.1f %%', refFrameG, Th(1) + Th(2), Th(3)));
        subplot(2,2,4);  imshow(labeloverlay(photoW, rustFaint & photoValid, 'Colormap', [1 0.6 0], 'Transparency', 0.55));  hold on;
        if any(hot(:)),  visboundaries(hot,  'Color', 'g', 'LineWidth', 0.8); end
        if any(bigM(:)), visboundaries(bigM, 'Color', 'c', 'LineWidth', 2); end
        if any(cold(:)), visboundaries(cold, 'Color', 'r', 'LineWidth', 0.8); end
        title('orange fill = rust colour (faint level), green = hot anomalies, cyan = largest hot region, red = cold');
        for ax = findall(f, 'Type', 'axes')', try, ax.Toolbar = []; catch, end; end
        stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        try
            exportgraphics(f, fullfile(outFolder, ['photo_validation_' pl.name '_' gt maskTag '_' stamp '.png']), 'Resolution', 200);
        catch ME
            fprintf('Figure for %s NOT saved (%s).\n', pl.name, ME.message);
        end
    end
end

%8.Results table
if ~isempty(summ)
    T = cell2table(summ, 'VariableNames', {'Plate','Group','Frames','RoiRustStrictPct','RoiRustLoosePct','RoiRustFaintPct','ThermHotPct','ThermColdPct','ThermLargestRegionPct','RhoLargeScale','RustDensityPct','PRustGivenAnomalyPct', ...
        'Lift','Rho','FramesWithBigRegion','BigRegionPctOfPlate','PrecisionLoose','RecallLoose','DiceLoose','PrecisionFaint','RecallFaint','DiceFaint'});
    writetable(T, fullfile(outFolder, ['photo_validation_results' maskTag '.csv']));
    fprintf('\nSaved: %s\n', fullfile(outFolder, ['photo_validation_results' maskTag '.csv']));
end

% ---------------------------------------------------------------------
function [rs, rl, rf] = rust_colour_lab(rgb, useMask, thr)
% Rust = Lab a* and b* both above the thresholds (strict / loose / faint).
lab = rgb2lab(im2double(rgb));
Af = imgaussfilt(lab(:,:,2), 1.2);
Bf = imgaussfilt(lab(:,:,3), 1.2);
rs = (Af > thr(1)) & (Bf > thr(2)) & useMask;
rl = (Af > thr(3)) & (Bf > thr(4)) & useMask;
rf = (Af > thr(5)) & (Bf > thr(6)) & useMask;
end

function [rs, rl, rf] = rust_colour_dark(rgb, useMask, offs, aMin)
% Rust = darker than the plate's bare steel (Otsu split of L* over the plate) and leaning red (a* > aMin).
% Loose / faint levels also take pixels clearly darker than the local steel (fine pits, thin streaks).
lab = rgb2lab(im2double(rgb));
Lf = imgaussfilt(lab(:,:,1), 0.8);
Af = imgaussfilt(lab(:,:,2), 0.8);
t = 100 * graythresh(Lf(useMask) / 100);
Ls = imresize(Lf, 0.25);
Lbg = imresize(medfilt2(Ls, [21 21], 'symmetric'), size(Lf), 'bilinear');      % local steel level
fineL = (Lf < Lbg - 6) & (Af > -1) & (Lf < t + 14);
fineF = (Lf < Lbg - 4) & (Af > -1) & (Lf < t + 14);
rs = (Lf < t + offs(1)) & (Af > aMin(1)) & useMask;
rl = (((Lf < t + offs(2)) & (Af > aMin(2))) | fineL) & useMask;
rf = (((Lf < t + offs(3)) & (Af > aMin(3))) | fineF | ((Lf < t + 10) & (Af > -1))) & useMask;
end

function bigM = big_regions(mask, rad, minA)
% Anomalous pixels in joined regions of at least minA pixels (same joining as thermal_anomaly_map).
cc = bwconncomp(imclose(mask, strel('disk', rad)));
a = cellfun(@numel, cc.PixelIdxList);
bigM = false(size(mask));
for i = find(a >= minA)
    bigM(cc.PixelIdxList{i}) = true;
end
bigM = bigM & mask;
end

function y = nan_smooth(x, valid, sg)
% Gaussian smoothing that ignores pixels outside 'valid'.
z = x;  z(~valid) = 0;
w = imgaussfilt(double(valid), sg);
y = imgaussfilt(z, sg) ./ max(w, eps);
y(~valid) = NaN;
end
