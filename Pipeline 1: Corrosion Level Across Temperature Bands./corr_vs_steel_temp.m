root = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\FINAL PIPELINES\corr pp';
addpath(root);
statsFolder = fullfile(root,'stats');
plates   = {'LIGHT','HALF'};
refName  = 'BARE NEW';  refExclude = [2146 2152 2153 2123];
bandLo   = [22 24 28 30 38 40 48 50];
maxFr = 12; sm = 2; pLow = 20;

fprintf('Reading bare reference...\n');
Rf = readFolder(fullfile(root,refName), sm, statsFolder);
Rf = Rf(~ismember([Rf.frame], refExclude));

rows = {};
for b = bandLo
    ir = pick(find([Rf.band]==b & [Rf.ok]), maxFr);
    if numel(ir) < 2, fprintf('band %d: <2 bare frames, skipped\n', b); continue; end
    rows(end+1,:) = {refName, b, numel(ir), median([Rf(ir).std]), NaN, NaN, NaN}; %#ok<SAGROW>
end
for p = 1:numel(plates)
    fprintf('Reading %s...\n', plates{p});
    P = readFolder(fullfile(root,plates{p}), sm, statsFolder);
    for b = bandLo
        ir = pick(find([Rf.band]==b & [Rf.ok]), maxFr);
        if numel(ir) < 2, continue; end
        DELTA = median([Rf(ir).spread]);
        it = pick(find([P.band]==b & [P.ok]), maxFr);
        if isempty(it), continue; end
        dT = nan(numel(it),1); hp = dT; sd = dT;
        for i = 1:numel(it)
            tS = P(it(i)).tS; mask = ~isnan(tS); v = tS(mask);
            hot = mask & (tS > prctile(v,pLow) + DELTA);
            hp(i) = 100*nnz(hot)/nnz(mask);
            if nnz(hot) > 0 && nnz(mask & ~hot) > 0
                dT(i) = mean(tS(hot)) - mean(tS(mask & ~hot));
            end
            sd(i) = P(it(i)).std;
        end
        f = median(hp,'omitnan')/100;
        model = sqrt(f*(1-f)) * median(dT,'omitnan');
        rows(end+1,:) = {plates{p}, b, numel(it), median(sd,'omitnan'), median(dT,'omitnan'), median(hp,'omitnan'), model}; %#ok<SAGROW>
    end
end
T = cell2table(rows,'VariableNames',{'Plate','Band','nFrames','PlateStdC','RustMinusSteelC','HotAreaPct','ModelStdC'});
disp(T);
writetable(T, fullfile(statsFolder,'rust_steel_dT.csv'));

function S = readFolder(folder, sm, statsFolder)
% Frames flagged by check_roi_fit (stats\frame_bands_<folder>.csv) are marked not ok.
[~, fname] = fileparts(char(regexprep(string(folder),'[\\/]+$','')));
csv = fullfile(statsFolder, ['frame_bands_' fname '.csv']);
bad = [];
if isfile(csv)
    C = readtable(csv,'TextType','string');
    bad = C.frame(~startsWith(strtrim(C.status),"OK"));
    fprintf('  %s: %d frame(s) flagged by check_roi_fit.\n', fname, numel(bad));
else
    fprintf('  %s: %s not found, ROI flags NOT applied. Run check_roi_fit first.\n', fname, csv);
end
files = dir(fullfile(folder,'IV_*.xlsx'));
S = struct('frame',{},'band',{},'ok',{},'std',{},'spread',{},'tS',{});
for k = 1:numel(files)
    t = regexp(files(k).name,'^IV_(\d+)\.xlsx$','tokens','once'); if isempty(t), continue; end
    fr = str2double(t{1});
    try
        d = load_thermal_matrix(fullfile(folder,files(k).name));
        [roi,info] = extract_plate_roi(d,false,'',true);
        mask = ~isnan(roi); z = roi; z(~mask) = median(roi(mask));
        tS = imgaussfilt(z,sm); tS(~mask) = NaN; v = tS(mask);
        n = numel(S)+1;
        S(n).frame  = fr;
        S(n).band   = floor(info.median/2)*2;
        S(n).ok     = nnz(mask)/numel(roi) >= 0.5 && ~ismember(fr, bad);
        S(n).std    = info.std;
        S(n).spread = prctile(v,95) - prctile(v,20);
        S(n).tS     = tS;
    catch
    end
end
end

function idx = pick(idx, nMax)
if numel(idx) > nMax, idx = idx(unique(round(linspace(1,numel(idx),nMax)))); end
end