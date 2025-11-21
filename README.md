# Multisine
This project contains code for generation of multisine signals and analysis of signals measured from tests that use multisines as excitation.

The multisine generation code is based on theory from "System Identification: A Frequency Domain Approach" by Pintelon and Schoukens. 

A multisine signal is a periodic signal. It consists of a mixture of sines of equal amplitude upto a desired bandwidth at the frequency resolution determined by the sampling frequency and the number of samples in one period. For example, if the sampling frequency is 1024Hz and the period is 4 seconds (= 4096 samples), then the frequency resolution = 1024Hz/4096 = 0.25Hz. The multisine of bandwidth, for example 100HZ, consists of sines of frequencies 0.25Hz, 0.5Hz, ..., 99.5Hz, 99.75Hz, 100Hz, i.e. 400 sine components all of equal amplitude. The phasess of these sines are adjusted such that the crest factor of the signal, i.e. the ratio of the highest peak to the lowest peak, is minimized to be close to 1.0. This ensures uniform amplitude of the signal over time and thus a continuously low signal noise ratio. This signal has been found to be well-suited from frequency response measurements. The periodicity of the input (and hence of the output) allows for application of FFT-based differentiation, integration and filtering during post-processing of data.

The multisine generation code is written in C and has a MEX wrapper, so it can be called from MATLAB.

To build the MEX function, call

```mex -D__MATLAB msingen.c fft.c```

The MATLAB function ```WindowedMsinGen``` is a higher-level wrapper that concatenates several periods of the multisine signal as well as includes ramp-up and ramp-down windows.

For analysis, data from windows where steady state has been achieved must be used.

The function ```multisineFRF``` can be used to perform frequency response analysis on measured data.
