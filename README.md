# Comparisons between SQUID- and OPM-based MEG systems
Research question: Would OPM-MEG systems in the long run be able to replace SQUID/Cryogenic/conventional MEG systems? If yes, how many OPM-MEG sensors do we need to outperform the existing 275-channel SQUID-MEG system from CTF?

My master thesis project is split in 3 parts: 1. simulations 2. phantom measurements 3. measuremts in human subjects

# Simulations

The code is based on the code provided by Robert Oostenveld in the course Neuroimaging II - Electrophysiological Methods and it is part of Konstantinos Tsilimparis' internship project.

The opm-meg system we try to simulate is the one ordered in the DCCN from FieldLine Inc. It has 32 sensors, 2 channels each (one radial, one tangential)

*forward_model.m*:

What I try to simulate: 
1. opm-meg closer to the scalp (like the eeg)
2. opm-meg have also 2 orthogonal channels per sensor (one radial, one tangential)

What I want to test:
1. Play with noise on the dipole (correlated vs uncorrelated, opm_sensornoise~2 *squid_sensornoise), 
2. Play with the number of opm sensors (32 sensors vs 64) & their location (on top of dipole vs homogeneous whole-head) 
3. Detection of 1 dipole vs 2  dipoles vs n dipoles
4. Superficial vs deep dipole

Metrics I want to use to quantify my comparsion:
1. Sensor-level: 
maximum field detected

2. Source-level:
a. dipole location error,
b. compare Vdata (ideal topography detected following dipole location) that can be explained by the Vmodel (measured topography following dipole location). For this, I can use

R-squared (coefficient of determination) - How good is the fit (measures the proportion of the variance in the Vdata that can be explained by the Vmodel.): R_squared = 1 - sum((Vdata_N20 - Vmodel_N20).^2) / sum((Vdata_N20 - mean(Vdata_N20)).^2)

Root Mean Square Error (RMSE): rmse = sqrt(mean((data - model).^2));
