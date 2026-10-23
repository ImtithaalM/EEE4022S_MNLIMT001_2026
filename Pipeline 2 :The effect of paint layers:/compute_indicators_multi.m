function out = compute_indicators_multi(baseROIs, targetROI, label, quiet)
% COMPUTE_INDICATORS_MULTI  Compares a target ROI against several bare ROIs and reports each indicator as mean and sd over the bare frames.
% Comparing against each bare frame separately (not an averaged image) keeps the spread ratio like-with-like; the sd is a built-in noise estimate.
% out = compute_indicators_multi(baseROIs, targetROI, label, quiet) returns out.mean.<field>, out.sd.<field>, out.fields, out.n, out.vals.

if nargin < 3 || isempty(label), label = 'target'; end
if nargin < 4 || isempty(quiet), quiet = false; end

%1.Indicators per bare frame
fields = {'offset','dT_mean','dT_max','dT_min','T_tilde', ...
          'N_tilde','N_tilde_med','std_ratio','iqr_ratio', ...
          'std_ratio_dt','iqr_ratio_dt'};
n = numel(baseROIs);
vals = zeros(n, numel(fields));
for i = 1:n
    r = compute_indicators(baseROIs{i}, targetROI, '', true);
    for k = 1:numel(fields)
        vals(i,k) = r.(fields{k});
    end
end

%2.Mean and sd
for k = 1:numel(fields)
    out.mean.(fields{k}) = mean(vals(:,k));
    out.sd.(fields{k})   = std(vals(:,k));
end
out.fields = fields;
out.n = n;
out.vals = vals;

%3.Printout
if quiet, return; end
fprintf('\n--- Indicators: %s vs %d bare frames (mean +/- sd) ---\n', label, n);
fprintf('Temp   : offset %+.2f +/- %.2f degC\n', out.mean.offset, out.sd.offset);
fprintf('MNTD   : dT_mean %+.4f +/- %.4f | dT_max %+.4f +/- %.4f | dT_min %+.4f +/- %.4f\n', ...
    out.mean.dT_mean, out.sd.dT_mean, out.mean.dT_max, out.sd.dT_max, out.mean.dT_min, out.sd.dT_min);
fprintf('NPC    : N_tilde %+.3f +/- %.3f\n', out.mean.N_tilde, out.sd.N_tilde);
fprintf('NPCmed : N_tilde_med %+.3f +/- %.3f\n', out.mean.N_tilde_med, out.sd.N_tilde_med);
fprintf('Spread : std ratio %.3f +/- %.3f | IQR ratio %.3f +/- %.3f\n', ...
    out.mean.std_ratio, out.sd.std_ratio, out.mean.iqr_ratio, out.sd.iqr_ratio);
fprintf('SprdDT : std ratio %.3f +/- %.3f | IQR ratio %.3f +/- %.3f\n', ...
    out.mean.std_ratio_dt, out.sd.std_ratio_dt, out.mean.iqr_ratio_dt, out.sd.iqr_ratio_dt);
end
