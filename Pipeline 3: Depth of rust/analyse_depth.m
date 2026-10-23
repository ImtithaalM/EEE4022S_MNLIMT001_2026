% ANALYSE_DEPTH  dT of each depth rectangle vs rust thickness (micrometer before - after), per temperature band, for the 2 mm and 4 mm plates.
% Reads stats\depth_frame_bands_<folder>.csv from check_depth_roi_fit (status OK, band not OUT). Active spots: 2 mm TM, BM; 4 mm TR2, B.
% Output: printed tables, two figures and stats\depth_vs_thickness_new.csv (analyse_depth_combined adds the old four spots).

clear plates res          % scripts share the workspace

%baseDir = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\rust pipeline';
baseDir = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\rust pipeline';

%1.Settings
% micrometer readings (mm), three per spot, rows in rectNames order; before = with rust, after = bare steel (NaN = left off the plots)
plates(1).name      = '2mm rust';
plates(1).csvName   = '2mm new';         % folder name used by check_depth_roi_fit
plates(1).marker    = 'o';
plates(1).rectNames = {'TM', 'BM'};
plates(1).before    = [2.61 2.60 2.62;       % TM  (with rust)
                       2.56 2.58 2.58];      % BM
plates(1).after     = [2.51 2.50 2.50;       % TM  (bare steel)
                       2.52 2.50 2.50];      % BM
plates(2).name      = '4mm rust';
plates(2).csvName   = '4mm new';
plates(2).marker    = 's';
plates(2).rectNames = {'TR2', 'B'};
plates(2).before    = [3.54 3.51 3.53;       % TR2
                       3.45 3.46 3.47];      % B
plates(2).after     = [3.41 3.42 3.40;       % TR2
                       3.41 3.39 3.42];      % B

% original four spots (now shiny, not analysed), readings kept for reference
% 2 mm  before: TL 2.70 2.74 2.86 | TR 2.50 2.49 2.51 | BR 2.59 2.61 2.62 | BL 2.56 2.54 2.53   (BR and BL swapped back)
%       after : TL 2.50 2.51 2.51 | TR 2.50 2.49 2.51 | BR 2.53 2.50 2.51 | BL 2.51 2.49 2.50
% 4 mm  before: TL 3.41 3.39 3.40 | TR 3.59 3.59 3.56 | BR 3.61 3.61 3.60 | BL 3.48 3.50 3.45
%       after : TL 3.40 3.39 3.41 | TR 3.39 3.43 3.40 | BR 3.42 3.41 3.41 | BL 3.41 3.39 3.41

nP = numel(plates);
res = cell(1, nP);

%2.Rust thickness and dT per band
allBands = strings(0, 1);
for p = 1:nP
    pl = plates(p);
    rn = pl.rectNames;
    nR = numel(rn);
    pl.rustMm = mean(pl.before, 2) - mean(pl.after, 2);
    pl.rustSE = sqrt(var(pl.before, 0, 2)/3 + var(pl.after, 0, 2)/3);   % standard error of the difference

    fprintf('\n[%s] rust thickness per rectangle (mean before - mean after):\n', pl.name);
    for k = 1:nR
        fprintf('  %-4s before %.3f  after %.3f  rust %.3f mm (+/- %.3f)\n', rn{k}, ...
            mean(pl.before(k,:)), mean(pl.after(k,:)), pl.rustMm(k), pl.rustSE(k));
    end
    if any(isnan(pl.rustMm))
        warning('analyse_depth:noMicrometer', '[%s] micrometer readings missing for: %s', ...
            pl.name, strjoin(rn(isnan(pl.rustMm)), ', '));
    end
    fprintf('  bare steel readings (after): %.3f to %.3f mm, spread %.3f mm\n', ...
        min(mean(pl.after, 2)), max(mean(pl.after, 2)), max(mean(pl.after, 2)) - min(mean(pl.after, 2)));

    T = readtable(fullfile(baseDir, 'stats', ['depth_frame_bands_' pl.csvName '.csv']), 'TextType', 'string');
    for k = 1:nR
        if ~ismember([rn{k} '_dT'], T.Properties.VariableNames)
            error('analyse_depth:missingColumn', ['Column %s_dT not found in the %s CSV. Re-run ' ...
                'check_depth_roi_fit with the updated depth_geometry first.'], rn{k}, pl.name);
        end
    end
    use = startsWith(T.status, "OK") & T.band ~= "OUT";
    pl.bands = unique(T.band(use), 'stable');
    nB = numel(pl.bands);
    pl.dTmean = nan(nB, nR);  pl.dTsd = nan(nB, nR);  pl.nUse = zeros(nB, 1);
    for b = 1:nB
        sel = use & (T.band == pl.bands(b));
        pl.nUse(b) = nnz(sel);
        for k = 1:nR
            x = T.([rn{k} '_dT'])(sel);
            pl.dTmean(b,k) = mean(x);  pl.dTsd(b,k) = std(x);
        end
    end
    allBands = union(allBands, pl.bands);
    res{p} = pl;

    fprintf('\n[%s] dT (degC), mean +/- sd over usable frames:\n', pl.name);
    fprintf('  %-8s %3s', 'Band', 'n');
    for k = 1:nR, fprintf('  %-17s', sprintf('%s (%.2fmm)', rn{k}, pl.rustMm(k))); end
    fprintf('\n');
    for b = 1:nB
        fprintf('  %-8s %3d', pl.bands(b), pl.nUse(b));
        for k = 1:nR, fprintf('  %6.2f +/-%-5.2f   ', pl.dTmean(b,k), pl.dTsd(b,k)); end
        fprintf('\n');
    end
end

[~, ordB] = sort(str2double(extractBefore(allBands, '-')));
allBands = allBands(ordB);
%3.Figure 1: one line per band, both plates (circles/solid = 2 mm, squares/dashed = 4 mm)
figure('Name', 'dT vs rust thickness (all bands)', 'Color', 'w'); hold on; grid on;
cols = lines(numel(allBands));
lineStyles = {'-', '--'};
for bi = 1:numel(allBands)
    for p = 1:nP
        pl = res{p};
        b = find(pl.bands == allBands(bi), 1);
        if isempty(b), continue; end
        [xs, ord] = sort(pl.rustMm);
        h = errorbar(xs, pl.dTmean(b, ord)', pl.dTsd(b, ord)', 'Color', cols(bi,:), ...
            'LineStyle', lineStyles{p}, 'Marker', pl.marker, 'LineWidth', 1.5, 'MarkerSize', 7);
        if p == 1
            h.DisplayName = sprintf('%s degC', allBands(bi));
        else
            h.HandleVisibility = 'off';           % one legend entry per band
        end
    end
end
for p = 1:nP                                      % label the spots on the hottest band
    pl = res{p};
    b = find(pl.bands == allBands(end), 1);
    if isempty(b), continue; end
    for k = 1:numel(pl.rectNames)
        text(pl.rustMm(k) + 0.004, pl.dTmean(b,k) + 0.8, sprintf('%s (%dmm)', pl.rectNames{k}, 2*p), 'FontSize', 8);
    end
end
plot(nan, nan, 'k-o',  'DisplayName', '2 mm plate');
plot(nan, nan, 'k--s', 'DisplayName', '4 mm plate');
xlabel('Rust thickness (mm) = mean before - mean after');
ylabel('dT = rectangle median - plate reference median (degC)');
title('Depth of rust: 2 mm and 4 mm plates');
legend('Location', 'southeast');

%4.Figure 2: one panel per band
figure('Name', 'dT vs rust thickness', 'Color', 'w');
for bi = 1:numel(allBands)
    subplot(2, 2, bi); hold on; grid on;
    for p = 1:nP
        pl = res{p};
        b = find(pl.bands == allBands(bi), 1);
        if isempty(b), continue; end
        errorbar(pl.rustMm, pl.dTmean(b,:)', pl.dTsd(b,:)', pl.marker, 'LineWidth', 1.3, ...
            'MarkerSize', 7, 'DisplayName', sprintf('%s (n=%d)', pl.name, pl.nUse(b)));
        for k = 1:numel(pl.rectNames)
            text(pl.rustMm(k) + 0.004, pl.dTmean(b,k), pl.rectNames{k}, 'FontSize', 8);
        end
    end
    xlabel('Rust thickness (mm)');  ylabel('dT (degC)');
    title(sprintf('%s degC band', allBands(bi)));
    legend('Location', 'southeast');
end
sgtitle('dT of each rectangle vs rust thickness (rectangle median - plate reference median)');

%5.Save
out = table();
for p = 1:nP
    pl = res{p};
    nR = numel(pl.rectNames);
    rows = table(repmat(string(pl.name), nR, 1), string(pl.rectNames(:)), pl.rustMm, pl.rustSE, ...
        'VariableNames', {'plate', 'rect', 'rust_mm', 'rust_se_mm'});
    for bi = 1:numel(allBands)
        b = find(pl.bands == allBands(bi), 1);
        v = nan(nR, 1);
        if ~isempty(b), v = pl.dTmean(b,:)'; end
        rows.(matlab.lang.makeValidName("dT_" + allBands(bi))) = v;
    end
    out = [out; rows]; %#ok<AGROW>
end
writetable(out, fullfile(baseDir, 'stats', 'depth_vs_thickness_new.csv'));
fprintf('\nSaved: %s\n', fullfile(baseDir, 'stats', 'depth_vs_thickness_new.csv'));