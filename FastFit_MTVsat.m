function FastFit_MTVsat(inputs)
% =========================================================================
% FastFit_MTVsat
%
% Standalone MTVsat fitting from pre-aligned quantitative MRI maps and a single MTw scans.
%
% This function estimates PDsat, MTVsat, and derived maps from an MT-weighted
% acquisition.
%
% More details can be found here (see the fast implementation):
% https://www.biorxiv.org/content/10.1101/2024.08.09.606771v2.full
%
% This pipeline requires the mrQ and vistasoft MATLAB toolboxes.
% Please ensure both are added to your MATLAB path before running the code.
% -------------------------------------------------------------------------
% REQUIRED INPUT (structure "inputs")
% -------------------------------------------------------------------------
%
% All images MUST:
%   - Be in the SAME voxel space (same resolution, orientation, size)
%   - Be co-registered
%   - Be NIfTI files readable by readFileNifti
%
% ---- Quantitative maps ----
% inputs.R1_path        : [string] Path to R1 map (units: 1/sec)
% inputs.B1_path        : [string] Path to B1 map (unitless scaling, ~1 = nominal FA)
% inputs.Gains_path     : [string] Receive field / gain map (unitless)
% inputs.PD_path        : [string] Proton density map (arbitrary scaling, mrQ-style)
%
% ---- MT acquisition ----
% inputs.MT_path        : [string] MT-weighted image (same space)
% inputs.TR             : [scalar] Repetition time (seconds)
% inputs.FA_deg         : [scalar] Flip angle (degrees)
%
% ---- Masks ----
% inputs.BM_path        : [string] Brain mask (binary)
% inputs.WFmask_path    : [string] CSF mask (for calibration)
%
% ---- Optional maps (for derived outputs) ----
% inputs.TV_path        : [string] Tissue volume map (0–1)
% inputs.TVfull_path    : [string] Full tissue volume (for calibration)
%
% ---- Output ----
% inputs.outdir         : [string] Output directory for NIfTI maps
%
% -------------------------------------------------------------------------
% OUTPUTS (written to disk)
% -------------------------------------------------------------------------
% PDsat.nii.gz
% MTVsat.nii.gz
% R1sat.nii.gz
% SWF.nii.gz
% relaxivity_map.nii.gz
% And the same outputs following direct effect calibration:
% MTVsat_calibrated.nii.gz
% R1sat_calibrated.nii.gz
% SWF_map_calibrated.nii.gz
% relaxivity_map_calibrated.nii.gz
% nonMTV_R1_calibrated.nii.gz

% -------------------------------------------------------------------------
% IMPORTANT ASSUMPTIONS
% -------------------------------------------------------------------------
% - R1 is in [1/sec]
% - TR is in [seconds]
% - FA is provided in degrees (converted internally to radians)
% - B1 is a multiplicative scaling of FA (NOT percent!)
% - Gains is consistent with PD scaling (mrQ convention)
% - All maps are perfectly aligned
%
% =========================================================================

%% =========================
% 0. Setup
% =========================
if ~exist(inputs.outdir, 'dir')
    mkdir(inputs.outdir);
end

FA = inputs.FA_deg .* pi ./ 180; % convert to radians
TR = inputs.TR;

%% =========================
% 1. Load maps
% =========================
disp('Loading maps...')

nii = readFileNifti(inputs.R1_path);   R1 = double(nii.data);   xform = nii.qto_xyz;
nii = readFileNifti(inputs.B1_path);   B1 = double(nii.data);
nii = readFileNifti(inputs.Gains_path);Gains = double(nii.data);
nii = readFileNifti(inputs.PD_path);   PD = double(nii.data);
nii = readFileNifti(inputs.MT_path);   Smt = double(nii.data);

nii = readFileNifti(inputs.BM_path);   BM = nii.data > 0;
nii = readFileNifti(inputs.WFmask_path); WFmask = nii.data > 0;
nii = readFileNifti(inputs.TV_path); TV = double(nii.data);



%% =========================
% 2. Fit PDsat
% =========================
if ~exist(fullfile(inputs.outdir,'PDsat.nii.gz'))
    disp('Fitting PDsat...')

    ind = find(BM);
    sol = zeros(size(ind));

    parfor ii = 1:length(ind)
        idx = ind(ii);

        if PD(idx) <= 0 || R1(idx) <= 0 || ~isfinite(Smt(idx))
            sol(ii) = 0;
            continue;
        end

        options = optimoptions('lsqnonlin','Display','none');

        PDsatfun = @(x) Smt(idx) - ...
            (x .* Gains(idx) .* sin(B1(idx).*FA) .* ...
            (1 - exp(-(R1(idx).*(PD(idx)./x))*TR))) ./ ...
            (1 - cos(B1(idx).*FA) .* exp(-(R1(idx).*(PD(idx)./x))*TR));

        x0 = max(PD(idx), 0.1);

        try
            sol(ii) = lsqnonlin(PDsatfun, x0, 0, 5, options);
        catch
            sol(ii) = 0;
        end
    end

    PDsat = zeros(size(PD));
    PDsat(ind) = sol;
    PDsat(~isfinite(PDsat)) = 0;

    dtiWriteNiftiWrapper(PDsat, xform, fullfile(inputs.outdir,'PDsat.nii.gz'));
else
    nii = readFileNifti(fullfile(inputs.outdir,'PDsat.nii.gz')); PDsat = double(nii.data);

end
%% =========================
% 3. MTVsat
% =========================
if ~exist(fullfile(inputs.outdir,'MTVsat.nii.gz'))

    disp('Computing MTVsat...')

    Ttheor = 4.300; % seconds

    PDcsf = Smt ./ (Gains .* ...
        ((1 - exp(-TR./Ttheor)) .* sin(B1.*FA) ./ ...
        (1 - exp(-TR./Ttheor).*cos(B1.*FA))));

    % Calibration (CSF peak)
    [csfValues, csfDensity] = ksdensity(PDcsf(WFmask));
    CalibrationVal = csfDensity(csfValues == max(csfValues));
    CalibrationVal = 1 ./ CalibrationVal;

    WFsat = PDsat .* CalibrationVal(1);
    MTVsat = 1 - WFsat;
    MTVsat(MTVsat < 0) = 0;

    dtiWriteNiftiWrapper(MTVsat, xform, fullfile(inputs.outdir,'MTVsat.nii.gz'));
else
    nii = readFileNifti(fullfile(inputs.outdir,'MTVsat.nii.gz')); MTVsat = double(nii.data);
end
%% =========================
% 4. Derived maps
% =========================
if ~exist(fullfile(inputs.outdir,'relaxivty_map.nii.gz'))

    disp('Computing derived maps...')

    R1sat = R1 .* PD ./ PDsat;
    R1sat(~isfinite(R1sat)) = 0;
    dtiWriteNiftiWrapper(R1sat, xform, fullfile(inputs.outdir,'R1sat.nii.gz'));

    SWF = (MTVsat - TV) ./ (1 - TV);
    SWF(MTVsat == 0 | TV == 0) = 0;
    SWF(SWF > 1 | SWF < 0) = 0;
    SWF(~BM) = 0;
    dtiWriteNiftiWrapper(SWF, xform, fullfile(inputs.outdir,'SWF.nii.gz'));

    relaxivity_map = (R1sat - R1) ./ (MTVsat - TV);
    relaxivity_map(MTVsat == 0 | TV == 0) = 0;
    relaxivity_map(~BM) = 0;
    dtiWriteNiftiWrapper(relaxivity_map, xform, fullfile(inputs.outdir,'relaxivty_map.nii.gz'));
end

%% =========================
% 5. Calibrate for the direct effect
% =========================
if ~exist(fullfile(inputs.outdir,'nonMTV_R1_calibrated.nii.gz'))

    CalibrationValSat=median(MTVsat(WFmask==1 & MTVsat>0)); % MTVsat in the CSF
    CalibrationValnonSat=median(TV(WFmask==1)); % MTV in the CSF
    calbval=CalibrationValSat-CalibrationValnonSat; % direct effect estimation

    MTVsat_calb=double(MTVsat)-calbval;
    MTVsat_calb(MTVsat_calb<0)=0;
    dtiWriteNiftiWrapper(MTVsat_calb, xform, fullfile(inputs.outdir,'MTVsat_calibrated.nii.gz'));

    R1satcalb=R1.*(1-TV)./(1-MTVsat_calb);
    R1satcalb(MTVsat_calb==0 | TV==0)=0;
    dtiWriteNiftiWrapper(R1satcalb, xform, fullfile(inputs.outdir,'R1sat_calibrated.nii.gz'));

    relaxivity_map=(R1satcalb-R1)./(MTVsat_calb-TV);
    relaxivity_map(BM==0)=0;
    relaxivity_map(MTVsat_calb==0 | TV==0)=0;
    dtiWriteNiftiWrapper(relaxivity_map, xform, fullfile(inputs.outdir,'relaxivty_map_calibrated.nii.gz'));

    SWF=(MTVsat_calb-TV)./(1-TV);
    SWF(MTVsat_calb==0 | TV==0)=0;
    SWF(SWF>1)=0;
    SWF(SWF<0)=0;
    SWF(BM==0)=0;
    dtiWriteNiftiWrapper(SWF, xform, fullfile(inputs.outdir,'SWF_calibrated.nii.gz'));

    nonMTV_R1=R1-relaxivity_map.*TV;
    dtiWriteNiftiWrapper(nonMTV_R1, xform, fullfile(inputs.outdir,'nonMTV_R1_calibrated.nii.gz'));

end
end