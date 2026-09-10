# Why the unoptimized crest factor is about sqrt(2 ln N)

`msingen.c` starts each run from a multisine with random phases, then runs an
iterative clipping-and-refit loop (see `msingen.c`, the `ITNO`-iteration loop)
to push the crest factor down. This note derives the theoretical crest factor
*before* that optimization runs, i.e. for a plain random-phase multisine, and
explains why the optimization loop is worth having.

## Setup

A multisine with `N` equal-amplitude harmonics and random phases is

    x(t) = sum_{k=1}^{N} a_k * cos(2*pi*f_k*t + phi_k)

with `phi_k` i.i.d. uniform on `[0, 2*pi)` and all `a_k` equal.

## Step 1: the amplitude distribution is approximately Gaussian

For a fixed time `t`, `x(t)` is a sum of `N` independent-ish random terms (one
per harmonic, each with an independent random phase). By the central limit
theorem, for large `N` this sum is approximately Gaussian, with variance equal
to the signal's mean-square value (RMS^2). This is exactly why
`getTimeFunction()` in `msingen.c` normalizes the time-domain signal by its
RMS before computing a crest factor: after normalization, the signal behaves
like a standard Gaussian process, `x(t) ~ N(0, 1)`.

## Step 2: crest factor is a question about extremes of that Gaussian process

The crest factor is `max|x(t)| / rms(x)`. With `x(t)` normalized to unit
variance, this is just asking: what is the expected maximum of a Gaussian
process sampled some effective number of `M` "independent" times?

For `M` i.i.d. `N(0, 1)` draws, the classical Gaussian extreme-value
asymptotic gives, to leading order,

    E[max] ~= sqrt(2 * ln(M))

(there are lower-order correction terms involving `ln ln M`, but
`sqrt(2 ln M)` is the dominant term as `M` grows).

## Step 3: what stands in for M

A multisine is a smooth, bandlimited signal, not `M` truly independent
samples. But Rice's formula (the expected number of local maxima per unit
time of a stationary Gaussian process is set by its bandwidth) says a signal
built from `N` harmonics has on the order of `N` effectively independent
peaks over one period. So `M ~= N` = `numFreq`, giving

    crest factor ~= sqrt(2 * ln(numFreq))

## Numerical check

For `numFreq = 400`: `sqrt(2 * ln 400) ~= sqrt(11.98) ~= 3.46`.

This matches the very first (unoptimized) iteration measured by
instrumenting the MATLAB mex build directly (`crxopt` at `i = 0`):

```
i = 0, crxopt = 3.63357
```

The iterative clipping-and-refit loop then reduces this toward a converged
crest factor around 1.6-2.0 over ~200 iterations. A lower crest factor means
more of a peak-limited actuator's/DAC's dynamic range is used by genuine
signal content rather than being reserved for rare peaks -- i.e. better
signal-to-noise ratio for the same peak amplitude limit. This is the
motivation, from Pintelon & Schoukens ("System Identification: A Frequency
Domain Approach", also cited in the top-level `README.md`), for using
crest-factor-optimized rather than plain random-phase multisines as
excitation signals.
