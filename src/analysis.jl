using Statistics: Statistics, mean, var


"""
Check that an `AccelerationTrace` RMS is consistent with Monte-Carlo statistics.

```julia
n_sigma = check_acceleration_trace(psd, trace; verbose=true)
```

returns the deviation of the oscillator-amplitude-based RMS from the
band-limited PSD RMS, expressed in units of the Monte-Carlo standard deviation
σ. The band limits and oscillator count are read directly from `trace`. Values
≤ 4 are expected ~99.994% of the time for correct generation; larger values
suggest a bug. Whether `target_rms` was given during generation does not
matter: the raw RMS is recovered from the stored oscillator amplitudes.

# Keyword Arguments
- `verbose`: if `true`, print expected/raw RMS and σ summary; default `false`

## Monte-Carlo convergence

With M = `n_oscillators` sinusoids with amplitudes Aₖ = √(2 S(fₖ) Δf/M) and
frequencies fₖ ~ Uniform[f_band_min, f_band_max] (Δf = f_band_max − f_band_min),
the raw squared RMS is

    raw_rms² = Σₖ Aₖ²/2 = (Δf/M) Σₖ S(fₖ)

— Δf times the sample mean of S(fₖ).  By the law of large numbers this
converges to

    E[raw_rms²] = Δf · E[S] = ∫_band S(f) df = expected_rms²

where `expected_rms = rms(psd; f_min=f_band_min, f_max=f_band_max)`.

Note: `raw_rms` here is computed from the oscillator amplitudes, not from the
time-domain waveform.  The time-domain RMS equals √(Σₖ Aₖ²/2) only when every
oscillator completes many cycles ((1/N) Σₙ cos²(ωₖtₙ + φₖ) → 1/2).  For
oscillators with fₖ < df = 1/(N·dt), this time average depends heavily on the
random phase and can deviate far from 1/2.  Using the amplitude-based formula
keeps `raw_rms` a pure MC quantity so that σ_rms below gives the correct z-score
regardless of window length.

The M PSD samples S(fₖ) = Aₖ²·M/(2·Δf) are recovered from the stored
oscillator amplitudes, and their sample variance estimates Var[S(f)] directly.
The variance of the sample mean then gives

    Var[raw_rms²] ≈ (Δf²/M) · var(S(fₖ))

and first-order error propagation (δ√x ≈ δx / 2√x) yields

    σ_rms ≈ Δf · √(var(S(fₖ)) / M) / (2 · expected_rms)
"""
function check_acceleration_trace(
    psd::PSDData,
    trace::AccelerationTrace;
    verbose::Bool = false,
)
    f_band_min = trace.f_band_min
    f_band_max = trace.f_band_max
    M = size(trace.oscillators, 2)
    Δf = f_band_max - f_band_min
    expected_rms = rms(psd; f_min = f_band_min, f_max = f_band_max)

    # Recover the pre-normalization PSD samples from the stored oscillator
    # amplitudes: Aₖ_raw = Aₖ_stored / η, and Aₖ_raw = √(2·S(fₖ)·Δf/M), so
    # S(fₖ) = Aₖ_raw² · M / (2·Δf).
    A_raw = trace.oscillators[1, :] ./ trace.eta
    S_samples = A_raw .^ 2 .* (M / (2 * Δf))

    # Compute raw_rms from the oscillator amplitudes, not from the time-domain
    # waveform.  The time-domain RMS equals √(Σₖ Aₖ²/2) only when every
    # oscillator completes many cycles, so that (1/N) Σₙ cos²(ωₖtₙ + φₖ) → 1/2.
    # Oscillators with fₖ below df = 1/(N·dt) complete less than one full cycle;
    # for those, the time average of cos² depends heavily on the random phase φₖ
    # and can deviate far from 1/2.  If the PSD is large at low frequencies, these
    # few oscillators carry most of the energy and dominate the time-domain RMS,
    # adding variance that is not captured by σ_rms below.  Using the amplitude-
    # based formula raw_rms² = Σₖ Aₖ²/2 = (Δf/M) Σₖ S(fₖ) = Δf · mean(S_samples)
    # keeps raw_rms a pure MC quantity (function of the sampled frequencies only),
    # so the σ_rms formula below gives the correct z-score regardless of window length.
    raw_rms = sqrt(Δf * mean(S_samples))

    # Var[raw_rms²] = (Δf²/M) · var(S(fₖ)); σ_rms via first-order propagation.
    σ_rms = Δf * sqrt(var(S_samples) / M) / (2 * expected_rms)
    n_sigma = abs(raw_rms - expected_rms) / σ_rms

    if verbose
        println("Check Monte-Carlo convergence of acceleration trace:")
        println(" expected rms (band-limited PSD): ", expected_rms, " m/s²")
        println(" raw rms (pre-normalization)     : ", raw_rms, " m/s²")
        println(" σ_rms (Monte Carlo)             : ", σ_rms, " m/s²")
        println(" deviation                       : ", n_sigma, " σ")
        println(" (≤ 4σ expected ~99.994% of the time)")
    end

    return n_sigma
end


"""
Check that the spectrum RMS is consistent with the source acceleration trace.

```julia
rms_err = check_acceleration_spectrum(trace, spectrum; verbose=true)
```

returns `|spectrum.rms − trace.rms| / trace.rms`, the relative difference
between the RMS computed independently from the FFT (√(∑|F(f)|²·df)) and the
time-domain RMS stored in `trace`. By the normalization convention this should
be zero up to floating-point precision (typically < 1e-14); larger values
indicate a bug in the spectrum normalization.

# Keyword Arguments
- `verbose`: if `true`, print both RMS values and the relative error; default `false`
"""
function check_acceleration_spectrum(
    trace::AccelerationTrace,
    spectrum::AccelerationSpectrum;
    verbose::Bool = false,
)
    rms_err = abs(spectrum.rms - trace.rms) / trace.rms

    if verbose
        println("Check spectrum normalization (Parseval):")
        println(" rms from spectrum √(∑|F|²·df): ", spectrum.rms, " m/s²")
        println(" rms from time-domain trace   : ", trace.rms, " m/s²")
        println(" relative error               : ", rms_err)
        println(" (should be ≈ 0 up to floating-point precision)")
    end

    return rms_err
end
