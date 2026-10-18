% TUNE_RUST  Teaches the rust detector from clicked examples on your photos, then shows four panels at different biases to pick one by eye.
% Per photo: click the plate corners, then 6-8 spots of clean steel, then 8-12 spots of rust (dark, orange, grey-brown, pits), Enter after each.
% Positive bias = stricter. Saves stats\rust_protos_<name>.mat, which validate_with_photo uses for that plate (set its 'bias' there).

%root     = 'C:\Users\Imtithaal Manuel\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\matlab pipeline';
root    = 'C:\Users\imtit\OneDrive - University of Cape Town\4th year\sem 2\EEE4022S\matlab pipeline';
photoDir = fullfile(root, 'DIGITAL');
photos   = {'HALF','10mm_half.jpg'};
biasSets = struct('HALF', [-3 -2 -1 0]);
maxSide  = 1000;   % photo long side (px)
for p = 1:size(photos,1)
    biases = biasSets.(photos{p,1});
    I = imread(fullfile(photoDir, photos{p,2}));
    s = maxSide / max(size(I,1:2));  if s < 1, I = imresize(I, s); end
    lab = rgb2lab(im2double(I));

    %1.Click examples
    figure; imshow(I); title([photos{p,1} ': click plate corners, Enter']); [x, y] = getpts; close;
    pm = poly2mask(x, y, size(I,1), size(I,2));
    figure; imshow(I); title([photos{p,1} ': click 6-8 spots of CLEAN STEEL (pale and grey patches, spread over the plate), Enter']); [sx, sy] = getpts; close;
    figure; imshow(I); title([photos{p,1} ': click 8-12 spots of CORROSION: dark rust, orange, the dull grey-brown patches and fine pits (top right, bottom, left), Enter']); [rx, ry] = getpts; close;
    steelP = grab(lab, sx, sy);  rustP = grab(lab, rx, ry);

    %2.Show biases
    f = figure('Name', ['Rust: ' photos{p,1}], 'Position', [30 30 1300 900]);
    for k = 1:4
        [rs, ~] = rust_colour_rel(I, pm, steelP, rustP, biases(k), biases(k));
        subplot(2,2,k); imshow(labeloverlay(I, rs & pm, 'Colormap', [1 0 1], 'Transparency', 0.45));
        title(sprintf('%s: bias %g -> %.1f %% of plate', photos{p,1}, biases(k), 100*nnz(rs & pm)/nnz(pm)));
    end

    %3.Save
    save(fullfile(root, 'stats', ['rust_protos_' photos{p,1} '.mat']), 'steelP', 'rustP', 'pm');
end
function P = grab(lab, x, y)
P = zeros(numel(x), 3);
for i = 1:numel(x)
    cx = round(x(i)); cy = round(y(i));
    w = lab(max(cy-1,1):min(cy+1,end), max(cx-1,1):min(cx+1,end), :);   % 3x3 so fine dots keep their colour
    P(i,:) = squeeze(median(median(w,1),2))';
end
end
