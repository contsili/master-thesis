[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=contsili/master-thesis)

# A system with 128 OPMs can outperform a system with 275 SQUIDs 
This repository contains the publicly available code of my master's thesis. 

The work was conducted at the [Donders Institute for Brain, Cognition and Behaviour](https://www.ru.nl/donders/) and was supervised by [Robert Oostenveld](https://github.com/robertoostenveld).

## Content of the repository
- [readme.md](./readme.md) This readme
- [thesis.pdf](./thesis.pdf) The thesis .pdf-file
- [analysis/](./analysis) contains the code that has been used to perform all analyses and create the figures presented in this thesis
  - [characterization/](./analysis/characterization) Code to process and analyze the characterization measurements
  - [miscellaneous/](./analysis/miscellaneous) Code to generate statistics of the data and plot the filter settings
  - [offline/](./analysis/offline) Code for the offline analysis
  - [online/](./analysis/online) Code for the online analysis
  - [preprocessing/](./analysis/preprocessing) Code that was used to manually preprocess the participant data

## Short description
We study the wearable MEG sensors known as optically pumped magnetometers (OPMs). In comparison to the superconducting quantum interference devices (SQUIDs), OPM systems promise higher sensitivity and better spatial resolution. As a result, OPM systems have gained significant interest from clinics and neuroscience research institutions, including the Donders Institute. Since OPMs can be purchased and placed as individual sensors, clinics and neuroscience research institutions should choose how many sensors to get. For example, at the Donders Institute, we currently have 32 OPM sensors, but plan to get more in the future. This raises the question: how many OPM sensors should we get to outperform our current 275-sensor CTF SQUID system? This master thesis suggests that a system with 128 OPMs can outperform for cortical brain areas our current 275-sensor CTF SQUID system.


**Keywords:** OPMs, SQUIDs, MEG, source reconstruction, computational modelling