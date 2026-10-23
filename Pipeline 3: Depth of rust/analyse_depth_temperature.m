% ANALYSE_DEPTH_TEMPERATURE  How dT of every spot (old and new) changes with plate temperature: per-frame table, fit dT = a + b*T per spot, three figures.
% Run order: check_depth_roi_fit -> analyse_depth -> analyse_depth_combined -> analyse_depth_temperature. Needs the four depth_frame_bands CSVs and depth_vs_thickness_combined.csv.
% Output: stats\depth_all_spots_frames.csv, stats\depth_spot_sensitivity.csv; figures: slope vs thickness, normalised dT vs thickness, dT vs temperature.

%1.Settings
%baseDir = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\rust pipeline';
baseDir = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\rust pipeline';
statsDir = fullfile(baseDir, 'stats');

Tamb       = 20;               % ambient temperature (degC) for the normalised dT
yLimSlope  = [-0.15 0.10];     % plot window for the sensitivity (degC per degC)
yLimNorm   = [-0.40 0.20];     % plot window for the normalised dT
yLimDT     = [-4 2];           % plot window for dT in figure 3 (degC)
minFrames  = 4;                % fewest usable frames needed to fit a spot

plateNames = ["2mm rust", "4mm rust"];
markers    = {'o', 's'};
plotCol    = [0.00 0.45 0.74;  % 2 mm
              0.85 0.33 0.10]; % 4 mm
bands      = ["22-26", "28-32", "38-42", "48-52"];

% per-frame CSVs: plate, set, file, spots
src = {"2mm rust", "old", 'depth_frame_bands_2mm rust.csv', ["TL", "TR", "BR", "BL"];
       "2mm rust", "new", 'depth_frame_bands_2mm new.csv',  ["TM", "BM"];
       "4mm rust", "old", 'depth_frame_bands_4mm rust.csv', ["TL", "TR", "BR", "BL"];
       "4mm rust", "new", 'depth_frame_bands_4mm new.csv',  ["TR2", "B"]};

%2.Rust thickness per spot (analyse_depth_combined)
thickFile = fullfile(statsDir, 'depth_vs_thickness_combined.csv');
if ~isfile(thickFile)
    error('analyse_depth_temperature:noThickness', 'Run analyse_depth_combined first - %s not found.', thickFile);
end
TH = readtable(thickFile, 'TextType', 'string');
thKey = TH.plate + " " + TH.rect;

%3.Merge the per-frame CSVs
rows = table();
for s = 1:size(src, 1)
    f = fullfile(statsDir, src{s,3});
    if ~isfile(f)
        error('analyse_depth_temperature:missing', 'Per-frame CSV not found: %s', f);
    end
    T = readtable(f, 'TextType', 'string');
    n = height(T);
    usable = startsWith(T.status, "OK") & T.band ~= "OUT";
    for sp = src{s,4}
        nm = char(sp);
        R = table(repmat(src{s,1}, n, 1), repmat(src{s,2}, n, 1), T.frame, T.band, T.refMedian, ...
                  repmat(string(nm), n, 1), T.([nm '_median']), T.([nm '_dT']), usable, ...
                  'VariableNames', {'plate', 'set', 'frame', 'band', 'refMedian', 'spot', 'spotMedian', 'dT', 'usable'});
        rows = [rows; R]; %#ok<AGROW>
    end
end
[tf, loc] = ismember(rows.plate + " " + rows.spot, thKey);
rows.rust_mm    = nan(height(rows), 1);   rows.rust_mm(tf)    = TH.rust_mm(loc(tf));
rows.rust_se_mm = nan(height(rows), 1);   rows.rust_se_mm(tf) = TH.rust_se_mm(loc(tf));
rows.dTnorm = rows.dT ./ (rows.refMedian - Tamb);
if any(~tf)
    warning('analyse_depth_temperature:noThickness', 'No rust thickness for: %s (left out of the plots).', ...
        strjoin(unique(rows.plate(~tf) + " " + rows.spot(~tf))', ', '));
end
writetable(rows, fullfile(statsDir, 'depth_all_spots_frames.csv'));
fprintf('\nSaved: %s  (%d rows, %d spots)\n', fullfile(statsDir, 'depth_all_spots_frames.csv'), ...
    height(rows), height(unique(rows(:, {'plate', 'spot'}), 'rows')));

%4.Fit dT = a + b * refMedian per spot (usable frames)
U = unique(rows(:, {'plate', 'spot'}), 'rows', 'stable');
nS = height(U);
slope = nan(nS,1);  slopeSE = nan(nS,1);  icpt = nan(nS,1);  nFr = zeros(nS,1);
rust = nan(nS,1);   rustSE = nan(nS,1);   isNew = false(nS,1);
for i = 1:nS
    inSpot = rows.plate == U.plate(i) & rows.spot == U.spot(i);
    sel = inSpot & rows.usable;
    nFr(i) = nnz(sel);
    k = find(inSpot, 1);
    rust(i) = rows.rust_mm(k);  rustSE(i) = rows.rust_se_mm(k);  isNew(i) = rows.set(k) == "new";
    if nFr(i) >= minFrames
        [p, S] = polyfit(rows.refMedian(sel), rows.dT(sel), 1);
        Rinv = inv(S.R);
        C = (Rinv * Rinv') * S.normr^2 / S.df;
        slope(i) = p(1);  icpt(i) = p(2);  slopeSE(i) = sqrt(C(1,1));
    end
end
Sum = table(U.plate, U.spot, categorical(isNew, [false true], {'old', 'new'}), rust, rustSE, nFr, slope, slopeSE, icpt, ...
    'VariableNames', {'plate', 'spot', 'set', 'rust_mm', 'rust_se_mm', 'n_frames', 'slope_dT_per_degC', 'slope_se', 'intercept'});
writetable(Sum, fullfile(statsDir, 'depth_spot_sensitivity.csv'));

fprintf('\nSensitivity of dT to plate temperature, per spot:\n');
fprintf('  %-9s %-5s %-4s %8s %6s %10s %8s\n', 'Plate', 'Spot', 'Set', 'rust mm', 'n', 'slope', '+/-');
for p = 1:2
    idx = find(U.plate == plateNames(p));
    [~, o] = sort(rust(idx));
    for i = idx(o)'
        fprintf('  %-9s %-5s %-4s %8.3f %6d %10.4f %8.4f\n', U.plate(i), U.spot(i), string(Sum.set(i)), rust(i), nFr(i), slope(i), slopeSE(i));
    end
end

%5.Figure 1: sensitivity vs rust thickness
figure('Name', 'Sensitivity vs rust thickness', 'Color', 'w'); hold on; grid on;
for p = 1:2
    for flag = [false true]
        idx = find(U.plate == plateNames(p) & isNew == flag & ~isnan(rust) & ~isnan(slope))';
        if isempty(idx), continue; end
        on = slope(idx) >= yLimSlope(1);
        if flag, fc = plotCol(p,:); else, fc = 'none'; end
        if any(on)
            ii = idx(on);
            errorbar(rust(ii), slope(ii), slopeSE(ii), slopeSE(ii), rustSE(ii), rustSE(ii), markers{p}, ...
                'Color', plotCol(p,:), 'MarkerFaceColor', fc, 'MarkerSize', 8, 'LineWidth', 1.3, 'HandleVisibility', 'off');
            for q = ii
                text(rust(q) + 0.004, slope(q) + 0.006, U.spot(q), 'FontSize', 8);
            end
        end
        for q = idx(~on)                                 % bare steel: below the plot window
            text(0.004, yLimSlope(1) + 0.012 + 0.012*(p-1), sprintf('%s %s: %.2f (off scale)', extractBefore(plateNames(p), ' '), U.spot(q), slope(q)), ...
                'FontSize', 8, 'Color', plotCol(p,:));
        end
    end
end
h(1) = plot(nan, nan, 'o', 'Color', plotCol(1,:), 'MarkerFaceColor', 'none',      'DisplayName', '2 mm old');
h(2) = plot(nan, nan, 'o', 'Color', plotCol(1,:), 'MarkerFaceColor', plotCol(1,:), 'DisplayName', '2 mm new');
h(3) = plot(nan, nan, 's', 'Color', plotCol(2,:), 'MarkerFaceColor', 'none',      'DisplayName', '4 mm old');
h(4) = plot(nan, nan, 's', 'Color', plotCol(2,:), 'MarkerFaceColor', plotCol(2,:), 'DisplayName', '4 mm new');
legend(h, 'Location', 'southeast');
yline(0, 'k:', 'HandleVisibility', 'off');
ylim(yLimSlope);  xlim([-0.01 0.30]);
xlabel('Rust thickness (mm)');  ylabel('Sensitivity: dT per degC of plate temperature');
title('How strongly dT grows with temperature, vs rust thickness (bars = standard errors)');

%6.Figure 2: dT / (plate temperature - ambient) vs rust thickness, one line per band
figure('Name', 'Normalised dT vs rust thickness', 'Color', 'w');
bc = lines(4);
for p = 1:2
    subplot(1, 2, p); hold on; grid on;
    idx = find(U.plate == plateNames(p) & ~isnan(rust))';
    for b = 1:4
        xs = nan(numel(idx), 1);  ys = nan(numel(idx), 1);
        for j = 1:numel(idx)
            sel = rows.plate == U.plate(idx(j)) & rows.spot == U.spot(idx(j)) & rows.usable & rows.band == bands(b);
            if any(sel), xs(j) = rust(idx(j));  ys(j) = mean(rows.dTnorm(sel)); end
        end
        ok = ~isnan(ys);
        [xo, o] = sort(xs(ok));  yo = ys(ok);  yo = yo(o);
        plot(xo, yo, ['-' markers{p}], 'Color', bc(b,:), 'LineWidth', 1.3, 'MarkerSize', 6, 'DisplayName', sprintf('%s degC', bands(b)));
    end
    yline(0, 'k:', 'HandleVisibility', 'off');
    ylim(yLimNorm);  xlim([-0.01 0.30]);
    xlabel('Rust thickness (mm)');  ylabel(sprintf('dT / (T_{plate} - %g degC)', Tamb));
    title(plateNames(p));  legend('Location', 'southeast');
end
sgtitle('Normalised dT vs rust thickness (bands on top of each other = effect scales with temperature above ambient)');

%7.Figure 3: dT vs plate temperature, coloured by thickness
figure('Name', 'dT vs plate temperature', 'Color', 'w');
cm = parula(256);  rMax = max(rust);
for p = 1:2
    subplot(1, 2, p); hold on; grid on;
    for i = find(U.plate == plateNames(p))'
        sel = rows.plate == U.plate(i) & rows.spot == U.spot(i) & rows.usable;
        if isnan(rust(i)) || ~any(sel) || isnan(slope(i)), continue; end
        c = cm(min(256, 1 + round(255 * rust(i) / rMax)), :);
        x = rows.refMedian(sel);
        scatter(x, rows.dT(sel), 12, c, 'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
        xx = [min(x) max(x)];
        plot(xx, icpt(i) + slope(i) * xx, '-', 'Color', c, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        text(max(x) + 0.4, icpt(i) + slope(i) * max(x), sprintf('%s %.2f', U.spot(i), rust(i)), 'FontSize', 8);
    end
    colormap(gca, parula);  caxis([0 rMax]);  cb = colorbar;  cb.Label.String = 'Rust thickness (mm)';
    ylim(yLimDT);  xlim([20 58]);
    xlabel('Plate reference temperature (degC)');  ylabel('dT (degC)');
    title(plateNames(p));
end
sgtitle('dT vs plate temperature, all spots (labels: spot, rust thickness in mm; bare steel is off the bottom)');

fprintf('\nSaved: %s\n', fullfile(statsDir, 'depth_spot_sensitivity.csv'));
