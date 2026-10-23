function data = load_thermal_matrix(filepath)
% LOAD_THERMAL_MATRIX  Reads one IRSoft-exported thermal xlsx (480 x 640, degC) and checks it is usable.
% Warns on a wrong size, NaN cells or a flat (failed) export.
% Input: path to one .xlsx file. Output: 480 x 640 matrix, data(row, col) = temperature of that pixel.

data = readmatrix(filepath);

%1.Checks
[rows, cols] = size(data);
if rows ~= 480 || cols ~= 640
    warning('Unexpected size for %s (expected 480x640)', filepath);
end
if any(isnan(data(:)))
    warning('%s contains NaNs - %d cells', filepath, sum(isnan(data(:))));
end
if all(data(:) == data(1))
    warning('%s is a flat/constant matrix - likely a bad export', filepath);
end
end
