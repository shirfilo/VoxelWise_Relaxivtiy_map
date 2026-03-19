%% ========================================================================
% Run MTVsat, R1sat, relaxivity and SWF fitting from precomputed MPM maps, and a single MTw scan. 
%
% This script prepares the input structure and runs:
%    FastFit_MTVsat(inputs)
%
% All inputs must be co-registered and in the same space.
%
% More details can be found here (see fast implementation): 
% https://www.biorxiv.org/content/10.1101/2024.08.09.606771v2.full
%
% This pipeline requires the mrQ and vistasoft MATLAB toolboxes. 
% Please ensure both are added to your MATLAB path before running the code.
% ========================================================================

%% =========================
% 0. Setup paths
% =========================
addpath(genpath('/mrQ-master'))
addpath(genpath('/vistasoft-master'))

%% =========================
% 1. Subject directories
% =========================
base_dir = '';

hMRI_dir = fullfile(base_dir, 'hMRI_MPM');
mrQ_dir  = fullfile(base_dir, 'mrQ');

%% =========================
% 2. Build input structure
% =========================
inputs = struct();

% -------------------------------------------------------------------------
% Quantitative maps (same space!)
% -------------------------------------------------------------------------
tmp=dir(fullfile(hMRI_dir, 'Results', ...
    'LcpcaDenoised_*_R1.nii'));
inputs.R1_path = fullfile(tmp.folder,tmp.name) ;

inputs.B1_path    = fullfile(hMRI_dir, 'Results', 'Supplementary', ...
    'B1_reg.nii'); % unitless (≈1 = nominal FA)

inputs.Gains_path = fullfile(mrQ_dir, 'OutPutFiles_1', 'BiasMap', ...
    'Gains.nii.gz'); % receive field

inputs.PD_path    = fullfile(mrQ_dir, 'SPGR_1', 'RefIm_Align_1_1_1', ...
    'PD.nii.gz'); % proton density (mrQ scaling)

% -------------------------------------------------------------------------
% MT acquisition
% -------------------------------------------------------------------------
tmp=dir(fullfile(hMRI_dir, 'MPMCalc', ...
    '*_MTw_WLS1fit_TEzero.nii'));
inputs.MT_path = fullfile(tmp.folder,tmp.name) ;

inputs.TR     = 24 / 1000; % seconds
inputs.FA_deg = 6;         % degrees

% -------------------------------------------------------------------------
% Masks
% -------------------------------------------------------------------------
inputs.BM_path = fullfile(mrQ_dir, 'SPGR_1', 'RefIm_Align_1_1_1', ...
    'brainMask.nii.gz');

inputs.WFmask_path = fullfile(mrQ_dir, 'SPGR_1', 'RefIm_Align_1_1_1', ...
    'csf_seg_T1.nii.gz');

% -------------------------------------------------------------------------
% Optional maps
% -------------------------------------------------------------------------
inputs.TV_path = fullfile(mrQ_dir, 'OutPutFiles_1', 'BrainMaps', ...
    'TV_map.nii.gz');


% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
inputs.outdir = fullfile(base_dir, 'fitMTVsat');

%% =========================
% 3. Run fitting
% =========================
FastFit_MTVsat(inputs);
