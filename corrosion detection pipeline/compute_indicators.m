function out = compute_indicators(baseROI, targetROI, label, quiet)
% COMPUTE_INDICATORS  MNTD, NPC and spread indicators of a target plate ROI vs a bare-steel baseline ROI (both from extract_plate_roi, NaN = masked).
% Outputs: dT_mean/max/min and T_tilde (MNTD), N_tilde and N_tilde_med (NPC, baseline mean / median threshold), std_ratio and iqr_ratio (also _dt = slope removed).
% out = compute_indicators(baseROI, targetROI, label, quiet); label names the printout, quiet = true prints nothing.

b = baseROI(~isnan(baseROI));
t = targetROI(~isnan(targetROI));

Tbar_nc = mean(b);   TM_nc = max(b);   Tm_nc = min(b);
Tbar    = mean(t);   TM    = max(t);   Tm    = min(t);

out.base_mean = Tbar_nc;  out.base_max = TM_nc;  out.base_min = Tm_nc;
out.tgt_mean  = Tbar;     out.tgt_max  = TM;     out.tgt_min  = Tm;

%1.MNTD
out.dT_mean = (Tbar - Tbar_nc) / Tbar_nc;
out.dT_max  = (TM   - TM_nc)   / TM_nc;
out.dT_min  = (Tm   - Tm_nc)   / Tm_nc;
out.T_tilde = max([out.dT_mean, out.dT_max, out.dT_min]);

%2.NPC (baseline mean threshold)
out.fracBase = nnz(b > Tbar_nc) / numel(b);
out.fracTgt  = nnz(t > Tbar_nc) / numel(t);
out.N_tilde  = (out.fracTgt - out.fracBase) / out.fracBase;

out.offset = median(t) - median(b);   % overall temperature gap (degC)

%3.Spread
out.std_base  = std(b);    out.std_tgt  = std(t);
out.std_ratio = out.std_tgt / out.std_base;
out.iqr_base  = local_iqr(b);  out.iqr_tgt = local_iqr(t);
out.iqr_ratio = out.iqr_tgt / out.iqr_base;

%4.Spread after removing the overall slope
bd = local_detrend(baseROI);  bd = bd(~isnan(bd));
td = local_detrend(targetROI); td = td(~isnan(td));
out.std_base_dt  = std(bd);   out.std_tgt_dt = std(td);
out.std_ratio_dt = out.std_tgt_dt / out.std_base_dt;
out.iqr_ratio_dt = local_iqr(td) / local_iqr(bd);

%5.NPC (baseline median threshold)
Tmed_nc = median(b);
out.fracBaseMed = nnz(b > Tmed_nc) / numel(b);
out.fracTgtMed  = nnz(t > Tmed_nc) / numel(t);
out.N_tilde_med = (out.fracTgtMed - out.fracBaseMed) / out.fracBaseMed;

%6.Printout
if nargin < 3 || isempty(label), label = 'target'; end
if nargin < 4 || isempty(quiet), quiet = false; end
if quiet, return; end
fprintf('\n--- Indicators: %s vs bare baseline ---\n', label);
fprintf('Temp   : base median %.2f | target median %.2f | offset %+.2f degC\n', ...
    median(b), median(t), out.offset);
fprintf('MNTD   : dT_mean %+.4f | dT_max %+.4f | dT_min %+.4f  ->  T_tilde %+.4f\n', ...
    out.dT_mean, out.dT_max, out.dT_min, out.T_tilde);
fprintf('NPC    : base %.3f | target %.3f  ->  N_tilde %+.4f\n', ...
    out.fracBase, out.fracTgt, out.N_tilde);
fprintf('NPCmed : base %.3f | target %.3f  ->  N_tilde_med %+.4f\n', ...
    out.fracBaseMed, out.fracTgtMed, out.N_tilde_med);
fprintf('Spread : std %.3f -> %.3f (ratio %.3f) | IQR %.3f -> %.3f (ratio %.3f)\n', ...
    out.std_base, out.std_tgt, out.std_ratio, out.iqr_base, out.iqr_tgt, out.iqr_ratio);
fprintf('SprdDT : std %.3f -> %.3f (ratio %.3f) | IQR ratio %.3f\n', ...
    out.std_base_dt, out.std_tgt_dt, out.std_ratio_dt, out.iqr_ratio_dt);
if abs(out.offset) > 0.3
    fprintf('NOTE   : frames differ by %+.2f degC.\n', out.offset);
end
end

function q = local_iqr(v)
% Interquartile range (no toolbox needed).
v = sort(v(:));
n = numel(v);
q = v(max(1, round(0.75*n))) - v(max(1, round(0.25*n)));
end

function res = local_detrend(roi)
% Subtracts a least-squares plane (a + b*x + c*y) fitted to the valid pixels.
[nr, nc] = size(roi);
[X, Y] = meshgrid(1:nc, 1:nr);
X = (X - nc/2) / nc;  Y = (Y - nr/2) / nr;
m = ~isnan(roi);
A = [ones(nnz(m), 1), X(m), Y(m)];
coef = A \ roi(m);
res = roi - (coef(1) + coef(2)*X + coef(3)*Y);
end
