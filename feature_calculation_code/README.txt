This directory of MATLAB scripts contains the code necessary to calculate the feature values for a collection of experiments.  SEE calculate_features.m.

The code requires that the AllMeasurements.mat resulting from MACKtrack (Adelaja et al. 2021) has already been generated and is associated with an experimental ID using the google spreadsheet "Scope Runs" (directory information will need to be updated in locations.mat).
OR list of csv files with time series measurements as provided at Mendeley Data https://doi.org/10.17632/gkxzb5hcmk.1. 

The code has the following MATLAB requirements:

'MATLAB'	'9.7'
'Signal Processing Toolbox'	'8.3'
'Image Processing Toolbox'	'11.0'
'Statistics and Machine Learning Toolbox'	'11.6'
'Wavelet Toolbox'	'5.3'
'Curve Fitting Toolbox'	'3.5.10'

The code has the following direct dependicies on MACKtrack files:

MACKtrack\CellQuery\subfunctions\alignTrajectories.m
MACKtrack\CellTrack\subfunctions\backgroundcalculate.m
MACKtrack\CellTrack\subfunctions\checkread.m
MACKtrack\Run\readScopeRuns.m
MACKtrack\Utilities\globalpeaks.m
MACKtrack\Utilities\maketicks.m
MACKtrack\Utilities\modebalance.m
MACKtrack\Utilities\namecheck.m
MACKtrack\Utilities\quickthresh.m
MACKtrack\Utilities\smoothrows.m

MACKtrack can be dowloaded at https://github.com/Adewunmi91/MACKtrack