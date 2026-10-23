%% noise_check.m  -  camera noise floor from back-to-back frames, flat regions only
% sigma_noise = std(difference of two near-identical frames) / sqrt(2)
% Only the flattest 30% of pixels are used, so small frame-to-frame shifts
% at edges cannot inflate the estimate.

clear; clc;

%% ---- settings ----
folder   = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\FINAL PIPELINES\paint pp\BARE PAINT';
maxGap   = 2;       % max difference in frame number (IV_xxxxx)
maxDMed  = 0.3;     % max difference in plate median (C)
minCorr  = 0.98;    % min correlation between the two frames
flatPct  = 30;      % use the flattest this-many percent of pixels

% ROI (469 x 333)
r1 = 59 + 28;   r2 = 434 - 15;
c1 = 56 + 35;   c2 = 579 - 20;

%% ---- list frames in numeric order ----
files = dir(fullfile(folder, 'IV_*.xlsx'));
num = nan(numel(files),1);
for k = 1:numel(files)
    t = regexp(files(k).name, 'IV_(\d+)', 'tokens', 'once');
    if ~isempty(t), num(k) = str2double(t{1}); end
end
[num, idx] = sort(num);
files = files(idx);
nF = numel(files);
fprintf('Found %d frames in %s\n', nF, folder);

[X, Y] = meshgrid(1:(c2-c1+1), 1:(r2-r1+1));
A = [X(:) Y(:) ones(numel(X),1)];

%% ---- scan neighbouring pairs ----
sigMad = []; sigStd = []; sceneStd = []; used = {};
prev = []; prevNum = NaN; prevMed = NaN;

for k = 1:nF
    T = readmatrix(fullfile(folder, files(k).name));
    if ~isequal(size(T), [480 640])
        warning('%s is %dx%d, expected 480x640', files(k).name, size(T,1), size(T,2));
        prev = []; continue
    end
    R = T(r1:r2, c1:c2);
    med = median(R(:), 'omitnan');

    if ~isempty(prev) && (num(k) - prevNum) <= maxGap && abs(med - prevMed) <= maxDMed
        ok = isfinite(R(:)) & isfinite(prev(:));
        cc = corr(R(ok), prev(ok));
        if cc >= minCorr
            [gm, ~] = imgradient(imgaussfilt(prev, 1));        % local gradient
            flat = gm < prctile(gm(:), flatPct);               % flattest pixels
            ok2 = ok & flat(:);
            d = R(:) - prev(:);
            p = A(ok2,:) \ d(ok2);                             % remove slow drift (plane)
            d = d(ok2) - A(ok2,:) * p;
            sigMad(end+1,1)   = 1.4826 * mad(d, 1) / sqrt(2);  %#ok<SAGROW>
            sigStd(end+1,1)   = std(d) / sqrt(2);              %#ok<SAGROW>
            sceneStd(end+1,1) = std(R(ok));                    %#ok<SAGROW>
            used{end+1,1}     = sprintf('%d-%d', prevNum, num(k)); %#ok<SAGROW>
        end
    end
    prev = R; prevNum = num(k); prevMed = med;
end

%% ---- report ----
fprintf('\nAccepted pairs: %d\n', numel(sigMad));
if isempty(sigMad)
    fprintf('No valid pairs. Loosen maxGap (e.g. 5), maxDMed (0.5) or minCorr (0.95).\n');
    return
end
disp(table(used, 1000*sigMad, 1000*sigStd, sceneStd, ...
    'VariableNames', {'pair','noise_MAD_mK','noise_std_mK','plate_std_C'}));

ns = median(sigMad);
sc = median(sceneStd);
fprintf('Noise (robust, MAD): median %.1f mK  (range %.1f - %.1f mK)\n', ...
    1000*ns, 1000*min(sigMad), 1000*max(sigMad));
fprintf('Datasheet NETD:      50 mK\n');
fprintf('Plate std (scene):   median %.2f C\n', sc);
fprintf('Noise share of variance: %.3f %%\n', 100 * ns^2 / sc^2);