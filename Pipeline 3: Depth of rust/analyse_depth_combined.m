% ANALYSE_DEPTH_COMBINED  dT vs rust thickness for all six spots per plate: the four old spots (typed in below) plus the new spots from analyse_depth.
% Run order: check_depth_roi_fit -> analyse_depth -> analyse_depth_combined. The old and new spots were measured in different sessions.
% Output: printed table, two figures (per band, per plate) and stats\depth_vs_thickness_combined.csv. Bare-steel spots read far below yLim and are labelled instead.

%1.Settings
%baseDir = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\rust pipeline';
baseDir = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\rust pipeline';

yLim = [-4 2];                               % plot window for dT (degC)
xLim = [-0.01 0.30];                         % plot window for rust thickness (mm)

bandVars = {'dT_22_26', 'dT_28_32', 'dT_38_42', 'dT_48_52'};
bandLbl  = {'22-26',    '28-32',    '38-42',    '48-52'};
plateNames = ["2mm rust", "4mm rust"];
markers  = {'o', 's'};
plotCol  = [0.00 0.45 0.74;                  % 2 mm
            0.85 0.33 0.10];                 % 4 mm

%2.Old four spots (first analyse_depth run)
oldPlate = ["2mm rust"; "2mm rust"; "2mm rust"; "2mm rust"; "4mm rust"; "4mm rust"; "4mm rust"; "4mm rust"];
oldRect  = ["TL"; "TR"; "BR"; "BL"; "TL"; "TR"; "BR"; "BL"];
% 2 mm BR and BL thicknesses corrected here (BR 0.093, BL 0.043)
oldRust  = [0.26; 0; 0.0933333333333333; 0.0433333333333333; 0; 0.173333333333333; 0.193333333333333; 0.0733333333333333];
oldSE    = [0.0481894409826698; 0.00816496580927709; 0.0124721912892465; 0.0105409255338946; ...
            0.00816496580927727; 0.0156347191994114; 0.00471404520791022; 0.015986105077709];
oldDT    = [  0.291953086853024    0.449343754695015    0.783666338239401    1.25117444992066;     % 2mm TL
             -5.74600505828858    -8.94878218724176   -17.4728660583496    -25.1039624214172;      % 2mm TR (bare)
             -0.313490867614725   -0.554939270019516  -0.775771004813071   -1.20212650299072;      % 2mm BR
             -0.1630859375        -0.657306157625635  -2.13424164908273    -3.22310113906862;      % 2mm BL
             -5.32272787094118   -10.6625806263515    -19.2621530805315    -27.5777271270752;      % 4mm TL (bare)
             -0.019326496124275   -0.0308611733572673 -0.297139304024841   -0.0318910598754847;    % 4mm TR
              0.00053443908691015 -0.0635583968389431  0.127320425851009    0.0252147674560547;    % 4mm BR
             -0.30714931488039    -0.440989358084519  -0.182919638497471   -0.606449508666991];    % 4mm BL
old = table(oldPlate, oldRect, oldRust, oldSE, oldDT(:,1), oldDT(:,2), oldDT(:,3), oldDT(:,4), ...
    'VariableNames', {'plate', 'rect', 'rust_mm', 'rust_se_mm', bandVars{:}});
old.set = repmat("old", height(old), 1);

%3.New spots (from analyse_depth)
newFile = fullfile(baseDir, 'stats', 'depth_vs_thickness_new.csv');
if ~isfile(newFile)
    error('analyse_depth_combined:noNew', 'Run analyse_depth first - %s not found.', newFile);
end
new = readtable(newFile, 'TextType', 'string');
new = new(:, {'plate', 'rect', 'rust_mm', 'rust_se_mm', bandVars{:}});
new.set = repmat("new", height(new), 1);

T = [old; new];
if any(isnan(T.rust_mm))
    warning('analyse_depth_combined:noMicrometer', 'Spots without a rust thickness are left out: %s', ...
        strjoin(T.rect(isnan(T.rust_mm))', ', '));
    T = T(~isnan(T.rust_mm), :);
end

% spots left out of the plots and the saved table
exclude = strings(0, 1);                     % e.g. ["2mm rust BL"]
key  = T.plate + " " + T.rect;
drop = ismember(key, exclude);
if any(drop)
    fprintf('\nLeft out on purpose: %s\n', strjoin(key(drop)', ', '));
    T = T(~drop, :);
end

%4.Table
fprintf('\nAll spots: rust thickness and dT (degC)\n');
fprintf('  %-9s %-5s %-4s %8s %8s', 'Plate', 'Spot', 'Set', 'rust mm', '+/-');
for bi = 1:4, fprintf(' %9s', bandLbl{bi}); end
fprintf('\n');
for p = 1:2
    sel = find(T.plate == plateNames(p));
    [~, o] = sort(T.rust_mm(sel));
    for r = sel(o)'
        fprintf('  %-9s %-5s %-4s %8.3f %8.3f', T.plate(r), T.rect(r), T.set(r), T.rust_mm(r), T.rust_se_mm(r));
        for bi = 1:4, fprintf(' %9.2f', T.(bandVars{bi})(r)); end
        fprintf('\n');
    end
end

%5.Figure 1: one panel per band (hollow = old spot, filled = new spot)
figure('Name', 'dT vs rust thickness - all spots', 'Color', 'w');
for bi = 1:4
    subplot(2, 2, bi); hold on; grid on;
    for p = 1:2
        for isNew = [false true]
            sel = T.plate == plateNames(p) & ((T.set == "new") == isNew) & ~isnan(T.(bandVars{bi}));
            if ~any(sel), continue; end
            x = T.rust_mm(sel);  y = T.(bandVars{bi})(sel);  ex = T.rust_se_mm(sel);  nm = T.rect(sel);
            on = y >= yLim(1);                                    % points that fit inside the plot window
            if any(on)
                if isNew, fc = plotCol(p,:); else, fc = 'none'; end
                errorbar(x(on), y(on), zeros(nnz(on),1), zeros(nnz(on),1), ex(on), ex(on), markers{p}, ...
                    'Color', plotCol(p,:), 'MarkerFaceColor', fc, 'MarkerSize', 7, 'LineWidth', 1.2, ...
                    'HandleVisibility', 'off');
                for q = find(on)'
                    text(x(q) + 0.004, y(q) + 0.12, nm(q), 'FontSize', 8);
                end
            end
            for q = find(~on)'                                    % bare-steel spots: off the bottom of the plot
                text(max(x(q), 0) + 0.004, yLim(1) + 0.25 + 0.3*(p-1), ...
                    sprintf('%s %s: %.1f degC (off scale)', extractBefore(plateNames(p), ' '), nm(q), y(q)), 'FontSize', 8, 'Color', plotCol(p,:));
            end
        end
    end
    if bi == 1                                                    % legend entries (dummy points, so every combination shows)
        h(1) = plot(nan, nan, 'o',  'Color', plotCol(1,:), 'MarkerFaceColor', 'none',      'DisplayName', '2 mm old');
        h(2) = plot(nan, nan, 'o',  'Color', plotCol(1,:), 'MarkerFaceColor', plotCol(1,:), 'DisplayName', '2 mm new');
        h(3) = plot(nan, nan, 's',  'Color', plotCol(2,:), 'MarkerFaceColor', 'none',      'DisplayName', '4 mm old');
        h(4) = plot(nan, nan, 's',  'Color', plotCol(2,:), 'MarkerFaceColor', plotCol(2,:), 'DisplayName', '4 mm new');
        legend(h, 'Location', 'northwest');
    end
    xlim(xLim);  ylim(yLim);
    xlabel('Rust thickness (mm)');  ylabel('dT (degC)');
    title(sprintf('%s degC band', bandLbl{bi}));
end
sgtitle('dT of every spot vs rust thickness (hollow = old spots, filled = new spots)');

%6.Figure 2: one panel per plate, one line per band
figure('Name', 'dT vs rust thickness per plate', 'Color', 'w');
bc = lines(4);
for p = 1:2
    subplot(1, 2, p); hold on; grid on;
    for bi = 1:4
        sel = T.plate == plateNames(p) & ~isnan(T.(bandVars{bi}));
        if ~any(sel), continue; end
        [xs, o] = sort(T.rust_mm(sel));
        ys = T.(bandVars{bi})(sel);  ys = ys(o);
        plot(xs, ys, ['-' markers{p}], 'Color', bc(bi,:), 'LineWidth', 1.3, 'MarkerSize', 6, ...
            'DisplayName', sprintf('%s degC', bandLbl{bi}));
    end
    bare = T.plate == plateNames(p) & T.rust_mm < 0.005;               % bare-steel spot(s): far below the window
    if any(bare)
        v = T{bare, bandVars};
        text(0.005, yLim(1) + 0.3, sprintf('bare steel (0 mm): %.0f to %.0f degC (off scale)', max(v(:)), min(v(:))), 'FontSize', 8);
    end
    xlim(xLim);  ylim(yLim);
    xlabel('Rust thickness (mm)');  ylabel('dT (degC)');
    title(plateNames(p));
    legend('Location', 'southeast');
end
sgtitle('dT vs rust thickness, all spots per plate');

%7.Save
outFile = fullfile(baseDir, 'stats', 'depth_vs_thickness_combined.csv');
writetable(T, outFile);
fprintf('\nSaved: %s\n', outFile);