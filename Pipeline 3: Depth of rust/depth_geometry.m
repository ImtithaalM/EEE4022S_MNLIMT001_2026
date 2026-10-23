function g = depth_geometry(plate)
% DEPTH_GEOMETRY  Fixed plate outline and depth rectangles for the depth-of-rust test, g = depth_geometry('2mm') or ('4mm').
% Hard-coded from IRSoft line profiles (2 mm: IV_02393 / IV_02838, 4 mm: IV_02489 / IV_02919); redo only if the plate or camera moves. Coordinates are 1-based [x y] = [col row].
% Outline = four lines [x1 y1 x2 y2]; rectangle box = [colLeft colRight rowTop rowBottom]; g.activeIdx = rectangles that get statistics (the old four went shiny).

if nargin < 1 || isempty(plate), plate = '2mm'; end
plate = lower(char(plate));

switch plate
    case '2mm'
        % tilted plate
        g.top    = [ 49.88  31.08  588.03   67.08];      % P1 (539 samples)
        g.right  = [585.00  66.00  585.00  469.00];      % P2 (404 samples, vertical)
        g.bottom = [ 36.00 457.95  584.02  466.04];      % P3 (549 samples)
        g.left   = [ 37.84 458.12   51.25   30.75];      % P4 (428 samples)
        g.rectNames = {'TL', 'TR', 'BR', 'BL', 'TM', 'BM'};
        g.rectBox   = [204 316 137 250;                  % TL  (old, shiny now)
                       461 526 122 233;                  % TR  (old, shiny now)
                       469 535 310 424;                  % BR  (old, shiny now)
                        63 181 339 428;                  % BL  (old, shiny now)
                       202 336  81 126;                  % TM  new, top middle
                       208 298 371 451];                 % BM  new, bottom middle
        g.rectActive = [false false false false true true];
        % bare-steel spot (bottom centre), kept out of the reference region
        g.refExclude = [340 480 285 480];

    case '4mm'
        % nearly upright plate (fitted on IV_02489)
        g.top    = [ 58.0   54.0   580.0   95.0];        % 523 samples
        g.right  = [580.0   96.0   580.0  462.0];        % 367 samples, vertical
        g.bottom = [ 59.0  461.0   581.0  461.0];        % 523 samples, horizontal
        g.left   = [ 61.0   54.0    61.0  458.0];        % 405 samples, vertical
        g.rectNames = {'TL', 'TR', 'BR', 'BL', 'TR2', 'B'};
        g.rectBox   = [135 175 157 248;                  % TL  (cold: corrosion removed down to bare steel; old, shiny now)
                       328 443 154 237;                  % TR  (old, shiny now)
                       344 446 303 396;                  % BR  (old, shiny now)
                       134 240 321 421;                  % BL  (old, shiny now)
                       480 564 124 249;                  % TR2 new, top right
                       377 496 412 452];                 % B   new, bottom (cut at row 452 to stay inside the plate ROI)
        g.rectActive = [false false false false true true];
        % cold patch around TL, kept out of the reference region
        g.refExclude = [78 225 100 318];

    otherwise
        error('depth_geometry: unknown plate "%s" (use ''2mm'' or ''4mm'').', plate);
end
g.plate = plate;

% rectangles that get statistics
g.activeIdx   = find(g.rectActive);
g.activeNames = g.rectNames(g.activeIdx);
g.activeBox   = g.rectBox(g.activeIdx, :);

%Plate ROI inset
g.inset = [12 12 12 12];     % px inside the outline [top right bottom left], clear of edge blur

g.rectInset = 3;      % px trimmed off each rectangle side before statistics

%Hot-speck mask (as extract_plate_roi; cold pixels are kept, they are the signal)
g.kSigma     = 5;     % hot if more than kSigma robust sigma above the median
g.dilatePix  = 2;     % grow each masked pixel by this many pixels
g.speckWin   = 31;    % window for the speck test (px)
g.speckFrac  = 0.10;  % speck if less than this share of the window is hot

%Reference region
g.refMargin  = 5;     % px kept clear around each rectangle
end