# NoiseTAI.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ARLQCI.github.io/NoiseTAI.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ARLQCI.github.io/NoiseTAI.jl/dev/)
[![CI](https://github.com/ARLQCI/NoiseTAI/actions/workflows/CI.yml/badge.svg)](https://github.com/ARLQCI/NoiseTAI/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/ARLQCI/NoiseTAI.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ARLQCI/NoiseTAI.jl)

Generate synthetic random vibration time series with a prescribed acceleration [Power Spectral Density (PSD)](https://en.wikipedia.org/wiki/Spectral_density), for use in [Trapped Atom Interferometer (TAI)](https://arxiv.org/abs/2303.01100) noise simulations.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/ARLQCI/NoiseTAI.jl")
```

## Usage

Load a target PSD from a two-column ASCII file (Hz, a²/Hz), generate a
waveform with the matching spectral shape, and write the result:

```julia
using NoiseTAI:
    PSDData,
    AccelerationTrace,
    AccelerationSpectrum,
    write_acceleration_trace,
    write_acceleration_spectrum
using Random: Xoshiro

psd = PSDData("acc_plane.dat")

tlist = collect(range(0; step=6.5e-6, length=2^15))

# tlist[end] == 0.2129855 # seconds

trace = AccelerationTrace(
	psd;
    tlist,
    f_band_min = 2.05,
    f_band_max = 1025.0,
    n_oscillators = 50000,
    rng = Xoshiro(9)
)

spectrum = AccelerationSpectrum(trace)

write_acceleration_trace(trace, "trace.dat")
write_acceleration_spectrum(spectrum, "spectrum.dat")
```

See the [documentation](https://ARLQCI.github.io/NoiseTAI.jl/stable/) for the
full API reference and input file format.
