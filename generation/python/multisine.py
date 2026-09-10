"""Python wrapper around the msingen.c multisine generator.

Calls into libmsingen.so (built from msingen.c + fft.c with -D__PYCLIB;
see build.sh) via ctypes, so the numerics match the MATLAB mex build
exactly. windowed_msingen() ports WindowedMsinGen.m's post-processing
(window concatenation + raised-cosine ramp up/down) to numpy.
"""

import ctypes
from pathlib import Path

import numpy as np

_LIB_PATH = Path(__file__).resolve().parent / "libmsingen.so"

if not _LIB_PATH.exists():
    raise FileNotFoundError(
        f"{_LIB_PATH} not found. Build it first by running build.sh in "
        f"{_LIB_PATH.parent}."
    )

_lib = ctypes.CDLL(str(_LIB_PATH))
_lib.msingen.argtypes = [
    ctypes.c_ulong,
    ctypes.c_ulong,
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
]
_lib.msingen.restype = ctypes.c_long


def msingen(num_fft_points: int, num_freq: int) -> np.ndarray:
    """Generate a crest-factor-optimized multisine signal.

    num_fft_points -- number of samples in the multisine period.
    num_freq -- number of frequency components (bandwidth = num_freq / sampleRate).

    Returns a length-num_fft_points float64 numpy array.
    """
    multisinex = np.empty(num_fft_points, dtype=np.float64)
    _lib.msingen(num_fft_points, num_freq, multisinex)
    return multisinex


def windowed_msingen(
    num_fft_points: int,
    num_freq: int,
    num_windows: int,
    num_ramp_samples: int,
    scale: float,
    sample_rate: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Python port of WindowedMsinGen.m.

    Generates one multisine period via msingen(), repeats it num_windows
    times, and applies a raised-cosine ramp up/down over num_ramp_samples
    samples at the start/end, then scales by `scale`.

    Returns (t, wmsin): time vector in seconds and the windowed signal.
    """
    msin = msingen(num_fft_points, num_freq)

    wmsin = np.tile(msin, num_windows)

    ramp = np.ones_like(wmsin)
    ramp_up = (1 - np.cos(np.arange(num_ramp_samples) / (num_ramp_samples - 1) * np.pi)) / 2
    ramp[:num_ramp_samples] = ramp_up
    ramp[-num_ramp_samples:] = ramp_up[::-1]

    wmsin = wmsin * ramp * scale
    t = np.arange(len(wmsin)) / sample_rate

    return t, wmsin
