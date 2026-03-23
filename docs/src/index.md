# NoiseTAI


```@eval
using Markdown
using Pkg

VERSION = Pkg.dependencies()[Base.UUID("803bf123-8e03-4514-aad8-54507f43f3b1")].version

github_badge = "[![Github](https://img.shields.io/badge/ARLQCI-NoiseTAI.jl-blue.svg?logo=github)](https://github.com/ARLQCI/NoiseTAI.jl)"

version_badge = "![v$VERSION](https://img.shields.io/badge/version-v$(replace("$VERSION", "-" => "--"))-green.svg)"

Markdown.parse("$github_badge $version_badge")
```

Documentation for [NoiseTAI](https://github.com/ARLQCI/NoiseTAI.jl).

## Overview

`NoiseTAI` generates synthetic random vibration time series with a prescribed
acceleration [Power Spectral Density
(PSD)](https://en.wikipedia.org/wiki/Spectral_density), for use in Trapped
Atom Interferometer (TAI) noise simulations.

The primary workflow is:

1. Load a measured or analytically defined target PSD from a two-column ASCII
   file via [`PSDData`](@ref).
2. Generate a matching random waveform by superimposing `M` random sinusoids
   with physically derived amplitudes Aₖ = √(2·S(fₖ)·Δf/M) via
   [`AccelerationTrace`](@ref), optionally rescaled to a desired target RMS
   acceleration.
5. Write time-trace and spectrum files for use in physics simulations via
   [`write_acceleration_trace`](@ref) and
   [`write_acceleration_spectrum`](@ref).
