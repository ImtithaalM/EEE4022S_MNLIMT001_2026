% REFINE_RUST  Adds colour examples where the rust map is wrong and redraws it straight away.
% Each round shows the photo with the current rust in magenta: click missed rust, Enter, then magenta spots that are really steel, Enter.
% Repeat; a round with no clicks stops and saves.
root     = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\matlab pipeline';
photoDir = fullfile(root, 'DIGITAL');
plate    = 'PARTIAL';  photoFile = '10mm_partial.jpg';
bias     = -1;      % bias picked in tune_rust
maxSide  = 1000;    % photo long side (px)

%1.Load
file = fullfile(root, 'stats', ['rust_protos_' plate '.mat']);
Pz = load(file);  steelP = Pz.steelP;  rustP = Pz.rustP;  pm = Pz.pm;
I = imread(fullfile(photoDir, photoFile));
s = maxSide / max(size(I,1:2));  if s < 1, I = imresize(I, s); end
lab = rgb2lab(im2double(I));

%2.Refine
for round = 1:10
    [rs, ~] = rust_colour_rel(I, pm, steelP, rustP, bias, bias);
    figure('Name', 'Current rust map', 'Position', [30 30 1100 800]);
    imshow(labeloverlay(I, rs & pm, 'Colormap', [1 0 1], 'Transparency', 0.45));
    title(sprintf('Round %d: rust = magenta (%.1f %% of plate). Click MISSED rust, Enter', round, 100*nnz(rs & pm)/nnz(pm)));
    [mx, my] = getpts;
    title('Now click spots painted magenta that are really STEEL, Enter');
    [fx, fy] = getpts;  close;
    if isempty(mx) && isempty(fx), break; end
    rustP  = [rustP;  grab(lab, mx, my)];
    steelP = [steelP; grab(lab, fx, fy)];
end

%3.Save
save(file, 'steelP', 'rustP', '-append');
fprintf('Saved %d steel and %d rust examples to %s\n', size(steelP,1), size(rustP,1), file);
function P = grab(lab, x, y)
P = zeros(numel(x), 3);
for i = 1:numel(x)
    cx = round(x(i)); cy = round(y(i));
    w = lab(max(cy-1,1):min(cy+1,end), max(cx-1,1):min(cx+1,end), :);
    P(i,:) = squeeze(median(median(w,1),2))';
end
end
