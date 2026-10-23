function data = load_thermal_matrix(filepath)
% LOAD_THERMAL_MATRIX  Reads one IRSoft-exported thermal xlsx (480 x 640 matrix of degC, no header) and sanity-checks it.
% Warns on wrong size, NaN cells, or a flat matrix (likely a bad export). Call it on every raw frame before any other step.
% data = load_thermal_matrix('IV_00115.xlsx') returns data(row, col) in degC.

    data = readmatrix(filepath);

    %1.Size check
    [rows, cols] = size(data);
    fprintf('%s: %d x %d\n', filepath, rows, cols);
    if rows ~= 480 || cols ~= 640
        warning('Unexpected size for %s (expected 480x640)', filepath);
    end

    %2.NaN check
    if any(isnan(data(:)))
        warning('%s contains NaNs — %d cells', filepath, sum(isnan(data(:))));
    end

    %3.Flat matrix check
    if all(data(:) == data(1))
        warning('%s is a flat/constant matrix — likely a bad export', filepath);
    end
end
