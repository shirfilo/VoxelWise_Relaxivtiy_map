# VoxelWise_Relaxivtiy_map

MATLAB code for estimating MTVsat, R1sat, relaxivity, and SWF from quantitative MRI maps and a single MT-weighted scan.

Code accompanying the manuscript:
https://www.biorxiv.org/content/10.1101/2024.08.09.606771v2.full

# VoxelWise_Relaxivtiy_map

This pipeline requires the mrQ and vistasoft MATLAB toolboxes.
Ensure both are installed and added to your MATLAB path.

Usage

Edit and run:

run_FastFit_MTVsat.m

All inputs must be co-registered and in the same space.

Outputs

The following maps are saved:

PDsat.nii.gz

MTVsat.nii.gz

R1sat.nii.gz

SWF.nii.gz

relaxivity_map.nii.gz

Calibrated outputs (direct effect correction):

MTVsat_calibrated.nii.gz

R1sat_calibrated.nii.gz

SWF_map_calibrated.nii.gz

relaxivity_map_calibrated.nii.gz

nonMTV_R1_calibrated.nii.gz

Notes

Inputs must be in the same space (no registration is performed).

R1 in 1/sec, TR in seconds, FA in degrees.

For more details on the fast implementation, see the manuscript above.
