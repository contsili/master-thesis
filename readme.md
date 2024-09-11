[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=contsili/master-thesis)

# A system with 128 OPMs can outperform a system with 275 SQUIDs 
This repository contains the publicly available code of my master thesis. 

The work was conducted at the [Donders Institute for Brain, Cognition and Behaviour](https://www.ru.nl/donders/) and was supervised by [Robert Oostenveld](https://github.com/robertoostenveld).

## Content of the repository
- [master-thesis.pdf](./master-thesis.pdf) The thesis .pdf-file
- [analysis/](./analysis) contains the code used to perform all analyses and create the figures of this thesis
  - [coregistration/](./analysis/coregistration) code to perform co-registration of OPM and SQUID data.
  - [dmu/](./analysis/dmu) code to calculate the dipole moment uncertainty (DMU) for the OPM and SQUID data. DMU was calculated based on the somatosensory experiment and simulations (simulation of somatosensory activity & whole brain simulations). 
  - [dpu/](./analysis/dpu) code to calculate the dipole position uncertainty (DPU) for the OPM and SQUID data. DPU was calculated based on the somatosensory experiment and simulations (simulation of somatosensory activity & whole brain simulations). 
  - [empty-room/](./analysis/empty-room) code to analyze the empty room recordings of OPM and SQUID data.
  - [preprocessing/](./analysis/preprocessing) code to preprocess the OPM and SQUID data.
  - [supplementary/](./analysis/supplementary) code to create the supplementary material figures of this thesis.

## Short description
  We study the wearable MEG sensors known as optically pumped magnetometers (OPMs). In comparison to the superconducting quantum interference devices (SQUIDs), OPM systems promise higher sensitivity and better spatial resolution. As a result, OPM systems have gained significant interest from clinics and neuroscience research institutions, including the Donders Institute. Since OPMs can be purchased and placed as individual sensors, clinics and neuroscience research institutions should choose how many sensors to get. For example, at the Donders Institute, we currently have 32 OPM sensors, but plan to get more in the future. This raises the question: how many OPM sensors should we get to outperform our current 275-sensor CTF SQUID system? This master thesis suggests that a system with 128 OPMs can outperform for cortical brain areas our current 275-sensor CTF SQUID system.


**Keywords:** OPMs, SQUIDs, MEG, source reconstruction, computational modelling
