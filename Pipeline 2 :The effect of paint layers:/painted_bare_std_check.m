%% painted_bare_std_check.m  -  plate std per 2 degC band for the coated bare plates
root = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\FINAL PIPELINES\paint pp';   % same folder as painted_assess.m
addpath(root);
statsFolder = fullfile(root, 'stats');
plates    = {'BARE PAINT', 'BARE PAINT2'};
bandWidth = 2;
bandsWanted = [22 24 28 30 38 40 48 50];     % same bands as the Results tables
minUsable = 0.5;

rows = {};
for p = 1:numel(plates)
    folder = fullfile(root, plates{p});
    files = dir(fullfile(folder, 'IV_*.xlsx'));
    % frames flagged by check_roi_fit
    bad = [];
    csv = fullfile(statsFolder, ['frame_bands_' plates{p} '.csv']);
    if isfile(csv)
        C = readtable(csv, 'TextType', 'string');
        bad = C.frame(~startsWith(strtrim(C.status), "OK"));
        fprintf('%s: %d frame(s) flagged by check_roi_fit\n', plates{p}, numel(bad));
    else
        fprintf('%s: %s not found, ROI flags NOT applied\n', plates{p}, csv);
    end
    for k = 1:numel(files)
        t = regexp(files(k).name, '^IV_(\d+)\.xlsx$', 'tokens', 'once');
        if isempty(t), continue; end
        fr = str2double(t{1});
        if ismember(fr, bad), continue; end
        d = load_thermal_matrix(fullfile(folder, files(k).name));
        [roi, info] = extract_plate_roi(d, false, '', true);
        if nnz(~isnan(roi)) / numel(roi) < minUsable, continue; end
        rows(end+1,:) = {plates{p}, fr, info.median, info.std}; %#ok<SAGROW>
    end
end
R = cell2table(rows, 'VariableNames', {'Plate','Frame','MedianC','StdC'});
R.Band = floor(R.MedianC / bandWidth) * bandWidth;

fprintf('\n==== Plate std (degC) by plate and 2 degC band ====\n');
fprintf('%-12s %-6s %-6s %-7s %-7s %-7s %-9s %-9s\n', 'Plate','Band','n','MinStd','MedStd','MaxStd','Min/Med','Max/Med');
for p = 1:numel(plates)
    for b = bandsWanted
        m = strcmp(R.Plate, plates{p}) & R.Band == b;
        if ~any(m), continue; end
        s = R.StdC(m);
        fprintf('%-12s %d-%d  %-6d %-7.2f %-7.2f %-7.2f %-9.2f %-9.2f\n', plates{p}, b, b+bandWidth, ...
            numel(s), min(s), median(s), max(s), min(s)/median(s), max(s)/median(s));
    end
end
writetable(R, fullfile(statsFolder, 'painted_bare_std_all_frames.csv'));