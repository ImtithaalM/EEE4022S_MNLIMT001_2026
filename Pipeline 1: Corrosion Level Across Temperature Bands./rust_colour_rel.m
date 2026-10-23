function [rs, rl, info] = rust_colour_rel(rgb, useMask, steelProto, rustProto, biasS, biasL)
% RUST_COLOUR_REL  Nearest-example Lab classifier: a pixel is rust when it is closer to a clicked rust example than to any steel example.
% steelProto / rustProto are N x 3 [L a b] colours; score = dist(nearest steel) - dist(nearest rust).
% rs = strict (score > biasS), rl = loose (score > biasL), info.score = the score map.

%1.Score
lab = rgb2lab(im2double(rgb));
L = imgaussfilt(lab(:,:,1), 0.6);  A = imgaussfilt(lab(:,:,2), 0.6);  B = imgaussfilt(lab(:,:,3), 0.6);
px = [L(:) A(:) B(:)];
ds = min_dist(px, steelProto);
dr = min_dist(px, rustProto);
score = reshape(ds - dr, size(L));
score(~useMask) = -Inf;

%2.Rust masks (median filter fills JPEG-texture gaps)
rs = bwareaopen(medfilt2(score > biasS, [3 3]), 2);
rl = bwareaopen(medfilt2(score > biasL, [3 3]), 2);
info.score = score;
end

function d = min_dist(px, P)
% Distance from each pixel to its nearest example colour.
d = inf(size(px,1), 1);
for i = 1:size(P,1)
    di = sqrt(sum((px - P(i,:)).^2, 2));
    d = min(d, di);
end
end
