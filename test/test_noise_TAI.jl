using Test
using Printf
using Random
using TestingUtilities: @Test
import IOCapture

using NoiseTAI:
    AccelerationSpectrum,
    AccelerationTrace,
    PSDData,
    Ran2State,
    check_acceleration_spectrum,
    check_acceleration_trace,
    main,
    resample_acceleration_trace,
    rms

const PSD_FILE = joinpath(@__DIR__, "acc_plane.dat")
const FIXTURES = joinpath(@__DIR__, "fixtures")


"""
    normalize_floats(s; sigdigits=6, zero_threshold=1e-10)

Normalize a program output string for robust comparison by rounding all
floating-point literals to `sigdigits` significant digits using the `%g`
format, making comparisons insensitive to last-digit rounding variation.
Values with absolute value below `zero_threshold` are mapped to `"0"`, so that
floating-point noise values like `2.2e-16` compare equal to an exact `0`.
"""
function normalize_floats(s; sigdigits = 6, zero_threshold = 1e-10)
    s = replace(
        s,
        r"-?\d+\.\d+(?:[eE][+-]?\d+)?" =>
            m -> begin
                x = parse(Float64, m)
                abs(x) < zero_threshold ? "0" : @sprintf("%.*g", sigdigits, x)
            end
    )
    return s
end


"""
    parse_columns(file)

Parse a data file, skipping `#` comment lines.
Returns an N×C `Matrix{Float64}` (N rows = data points, C file columns).
"""
function parse_columns(file)
    rows = Vector{Float64}[]
    open(file) do fh
        for line in eachline(fh)
            s = strip(line)
            (isempty(s) || startswith(s, '#')) && continue
            push!(rows, parse.(Float64, split(s)))
        end
    end
    result = zeros(Float64, length(rows), length(rows[1]))
    for (i, row) in enumerate(rows)
        result[i, :] = row
    end
    return result
end


@testset "fixed seed small run" begin

    mktempdir() do tmpdir

        trace_file = joinpath(tmpdir, "trace.dat")
        spectrum_file = joinpath(tmpdir, "spectrum.dat")
        seed = -76989321
        # seed = -76989322
        n_oscillators = 50
        target_rms = 2.01956265143991
        # target_rms is arms_sample for acc_plane.dat with acc_unit=1
        log2_n_samples = 10
        n_samples = 2^log2_n_samples
        input_str = """
            log2 of n_samples : $log2_n_samples
            tstep (s): 1e-4
            PSD file: $PSD_FILE
            acc_unit: 1
            fpmin fpmax (Hz): 20 1025
            n_oscillators: $n_oscillators
            seed: $seed
            target_rms: $target_rms
            time trace file: $trace_file
            spectrum file: $spectrum_file
            resample: 2
            """
        captured = IOCapture.capture() do
            input_file = joinpath(tmpdir, "input.txt")
            write(input_file, input_str)
            open(input_file) do io
                redirect_stdin(io) do
                    main()
                end
            end
        end
        stdout_file = joinpath(tmpdir, "std.out")
        normalized_output = normalize_floats(captured.output)
        open(stdout_file, "w") do io
            write(io, normalized_output)
        end
        expected_output = normalize_floats("""
                                           # noise sample points in time = 2^m. Enter m<=20!
                                           log2_n_samples = $log2_n_samples
                                           n_samples = $n_samples
                                           tstep (s)?
                                           dt = 0.0001
                                           total_time = 0.1024
                                           f_max = 10000.0
                                           df = 9.765625
                                           Acceleration PSD sample file?
                                            (MUST BE sorted in ascending frequency)
                                           psd_file = "$PSD_FILE"
                                           Acceleration Unit G (m/s^2) = 9.81
                                           What is the unit for a, in m/s^2, in your sample?
                                           acceleration_unit = 1.0
                                           PSDData S(f) with 501 samples
                                              f_min: 2.05034 Hz
                                              f_max: 1025.169979 Hz
                                              rms:  $target_rms m/s²
                                           Noise band window (fpmin, fpmax)[Hz]?
                                             (Use something inside the above f_min, f_max)
                                           f_band_min = 20.0
                                           f_band_max = 1025.0
                                           Target number of oscillators nsmax?
                                                  Do < ntmax = $n_samples
                                           n_oscillators = $n_oscillators
                                            dt (s)             : 0.0001
                                            n_samples          : $n_samples
                                            total_time (s)     : 0.1024
                                            f_max (Hz)         : 10000.0
                                            Nyquist BW=0.5fmax : 5000.0
                                            df (Hz)            : 9.765625
                                            #Sample points     : 501
                                            f_min(sample, Hz)  : 2.05034
                                            f_max(sample, Hz)  : 1025.169979
                                            Noise band min (Hz): 20.0
                                            Noise band max (Hz): 1025.0
                                            # oscillators      : $n_oscillators

                                           Random seed? (int*4) < 0
                                           rng_seed = $seed
                                           Target RMS acceleration, accrms0, in m/s^2 ?
                                              Use the above reported arms_sample [m/s^2],
                                              or scale it up or down as you wish.
                                           target_rms = $target_rms
                                           Output file for time traces?
                                           time_trace_file = "$trace_file"
                                           Output file for spectrum?
                                           spectrum_file = "$spectrum_file"
                                           Initialize RAN2

                                           AccelerationTrace with $n_oscillators oscillators
                                              n_samples       : $n_samples
                                              dt (s)          : 0.0001
                                              f_band_min (Hz) : 20.0
                                              f_band_max (Hz) : 1025.0
                                              n_candidates    : 51
                                              rms (m/s²)      : $target_rms
                                              eta (target RMS / actual) : 1.5123091859600124

                                           Check Monte-Carlo convergence of acceleration trace:
                                            expected rms (band-limited PSD): 1.2762241779891113 m/s²
                                            raw rms (pre-normalization)     : 1.37083 m/s²
                                            σ_rms (Monte Carlo)             : 0.13857403962867668 m/s²
                                            deviation                       : 0.682676 σ
                                            (≤ 4σ expected ~99.994% of the time)

                                           AccelerationSpectrum with $n_samples samples
                                              df (Hz)  : 9.765625
                                              f_min (Hz): -5000.0
                                              f_max (Hz): 5000.0
                                              rms (m/s²): $target_rms

                                           Check spectrum normalization (Parseval):
                                            rms from spectrum √(∑|F|²·df): $target_rms m/s²
                                            rms from time-domain trace   : $target_rms m/s²
                                            relative error               : 2.198937525079677e-16
                                            (should be ≈ 0 up to floating-point precision)

                                           Output limited time sample? (1=Yes, else=no)
                                           do_resample = 2
                                           """)
        stdout_expected_file = joinpath(tmpdir, "std.expected.out")
        open(stdout_expected_file, "w") do io
            write(io, expected_output)
        end

        @test normalized_output == expected_output

        # --- output files exist ---
        @test isfile(trace_file)
        @test isfile(spectrum_file)

        # --- correct number of lines ---
        #   Header: 4 lines
        #     "# Acceleration time trace"
        #     "# Column 1: time (s)"
        #     "# Column 2: acceleration (m/s^2)"
        #     "# Column 3: displacement (m)"
        #   Data: one line per time step in `1:n_samples`
        @test countlines(trace_file) == 4 + n_samples

        # spectrum file:
        #   Header: 6 lines
        #     "# Two-sided acceleration spectrum"
        #     "# Normalized so that integral of PSD over all f equals RMS^2"
        #     "# Column 1: frequency (Hz)"
        #     "# Column 2: real part of spectral amplitude (m/s^2/sqrt(Hz))"
        #     "# Column 3: imaginary part of spectral amplitude (m/s^2/sqrt(Hz))"
        #     "# Column 4: two-sided power spectral density (m/s^2)^2/Hz"
        #   Data: one line per element of `1:n_freq_out`,
        #     n_freq_out = n_samples + 1.
        #     The +1 comes from the two-sided spectrum covering both positive
        #     and negative frequencies symmetrically: indices k = -n_half..n_half
        #     inclusive, giving n_samples/2 + 1 + n_samples/2 = n_samples + 1 bins.
        @test countlines(spectrum_file) == 6 + n_samples + 1

        trace = parse_columns(trace_file)
        spectrum = parse_columns(spectrum_file)

        a = trace[:, 2]   # acceleration (m/s²)
        d = trace[:, 3]   # displacement  (m)

        @testset "RMS normalization" begin
            # After generating the oscillator superposition, the program rescales
            # both waveforms by rms_scale = target_rms / acc_rms,
            # so that sqrt(sum(a²)/N) equals target_rms exactly.
            # We reproduce that formula here and check against target_rms.
            # The tolerance is 1e-8 rather than machine precision because the
            # trace file is written with "%15.8E" (9 significant digits),
            # and re-parsing introduces a relative rounding error of ~1e-9 per
            # sample that accumulates when squaring and summing.
            rms_val = sqrt(sum(a .^ 2) / length(a))
            @test rms_val ≈ target_rms rtol = 1e-8
        end

        @testset "displacement–acceleration consistency" begin
            # In the oscillator accumulation loop, each
            # sinusoid i contributes +acc_amplitude_i * cos(...) to a(t) and
            # -acc_amplitude_i/ω_i² * cos(...) to d(t). Both waveforms are then
            # scaled by the same rms_scale, so the relationship is preserved
            # after normalization.
            #
            # Denoting the i-th oscillator's contribution to a as aᵢ(t), we
            # have d(t) = -Σᵢ aᵢ(t)/ωᵢ², so the dot product sum(a .* d) = -Σᵢ
            # Σⱼ (1/ωⱼ²) sum(aᵢ .* aⱼ). The diagonal terms (i=j) are always
            # strictly negative; for approximately orthogonal sinusoids at
            # distinct frequencies the off-diagonal terms are small, so the
            # total is negative.
            @test sum(a .* d) < 0

            # Under the approximation that oscillator contributions are
            # orthogonal (valid for incommensurate frequencies over many
            # periods), rms(d)/rms(a) = sqrt(Σ Aᵢ²/ωᵢ⁴) / sqrt(Σ Aᵢ²),
            # which is bounded by the min-max inequality:
            #   1/ω_max² ≤ rms(d)/rms(a) ≤ 1/ω_min²
            # where ω_min = 2π·fpmin and ω_max = 2π·fpmax.
            # fpmin=20 Hz and fpmax=1025 Hz are the noise band limits in the input.
            f_min, f_max = 20.0, 1025.0
            ratio = sqrt(sum(d .^ 2) / sum(a .^ 2))
            @test ratio > 1 / (2π * f_max)^2
            @test ratio < 1 / (2π * f_min)^2
        end

        @testset "golden master" begin
            # The computation is fully deterministic given the fixed seed, so
            # output should be bit-for-bit identical across runs on the same
            # platform.  On the first run the fixtures are created from the
            # current output; thereafter every run is compared against them.
            # Uses a separate subdirectory from test_noise_TAI_orig.jl to
            # avoid collisions when both test files share the same FIXTURES path.
            fixtures = joinpath(FIXTURES, "julia")
            mkpath(fixtures)
            ref_trace = joinpath(fixtures, "trace.dat")
            if !isfile(ref_trace)
                cp(trace_file, ref_trace)
                @warn "Trace fixture not found; created from this run. Review and commit it."
            else
                @test trace == parse_columns(ref_trace)
            end
            ref_spectrum = joinpath(fixtures, "spectrum.dat")
            if !isfile(ref_spectrum)
                cp(spectrum_file, ref_spectrum)
                @warn "Spectrum fixture not found; created from this run. Review and commit it."
            else
                @test spectrum == parse_columns(ref_spectrum)
            end
        end

    end

end


@testset "API-based workflow with default RNG" begin

    # Parameters matching the interactive test, but driven entirely through
    # the Julia API.  A seeded Xoshiro instance (Julia's default RNG type)
    # makes the run deterministic without relying on the legacy Ran2State.
    dt            = 1e-4
    n_samples     = 2^10
    f_band_min    = 20.0
    f_band_max    = 1025.0
    n_oscillators = 50
    rng           = Random.Xoshiro(42)

    psd = PSDData(PSD_FILE; acceleration_unit = 1.0)

    tlist = collect(range(0.0, step = dt, length = n_samples))
    trace = AccelerationTrace(
        psd;
        tlist         = tlist,
        f_band_min    = f_band_min,
        f_band_max    = f_band_max,
        n_oscillators = n_oscillators,
        rng           = rng,
    )

    @testset "trace η is 1 (no rescaling)" begin
        # With target_rms=nothing the waveform is never rescaled, so η must be
        # exactly 1.0.
        @test trace.eta == 1.0
    end

    @testset "check_acceleration_trace within 4σ" begin
        # With target_rms=nothing, trace.rms is the raw Monte Carlo RMS.
        # check_acceleration_trace compares it to the band-limited PSD RMS and
        # expresses the deviation in units of the estimated MC standard deviation.
        # n_sigma ≤ 4 is expected ~99.994% of the time for any correct seed.
        n_sigma = check_acceleration_trace(psd, trace)
        @test n_sigma ≤ 4.0
    end

    spectrum = AccelerationSpectrum(trace)

    @testset "spectrum rms consistent with trace rms" begin
        # check_acceleration_spectrum returns the relative error between the
        # FFT-derived RMS and the time-domain RMS; Parseval's theorem guarantees
        # it is near zero up to floating-point precision (~1e-14).
        rms_err = check_acceleration_spectrum(trace, spectrum)
        @test rms_err < 1e-10
    end

    @testset "spectrum frequency axis" begin
        # The two-sided spectrum runs from -fmax/2 to +fmax/2 with spacing df.
        f_max = 1.0 / dt
        df = f_max / n_samples
        @test spectrum.freq_axis[1] ≈ -f_max / 2 rtol = 1e-12
        @test spectrum.freq_axis[end] ≈ f_max / 2 rtol = 1e-12
        @test length(spectrum.freq_axis) == n_samples + 1
        # Uniform spacing
        dfs = diff(spectrum.freq_axis)
        @test all(≈(df; rtol = 1e-12), dfs)
    end

    @testset "resample_acceleration_trace" begin
        # Resampling at a coarser grid preserves the oscillator parameters and
        # should give an RMS close to the original (within Monte-Carlo scatter).
        # With n_oscillators=50 the scatter is O(14%), so we use a loose bound.
        dt_coarse = 1e-3
        tlist_coarse = collect(range(0.0, step = dt_coarse, length = 512))
        resampled = resample_acceleration_trace(trace; tlist = tlist_coarse)
        @test resampled.oscillators === trace.oscillators
        # The resampled RMS should be within 50% of the original — generous
        # enough to be robust yet tight enough to catch total breakage.
        @test resampled.rms ≈ trace.rms rtol = 0.5
    end


end


@testset "rms function" begin

    psd = PSDData(PSD_FILE; acceleration_unit = 1.0)

    @testset "PSDData: no windowing matches psd.rms" begin
        # rms() with default f_min=0, f_max=Inf should reproduce the trapezoid
        # integral stored in psd.rms (the two formulations are mathematically
        # equivalent; any difference is pure floating-point rounding).
        @test rms(psd) ≈ psd.rms rtol = 1e-10
        @test rms(psd; f_min = 0.0, f_max = Inf) ≈ psd.rms rtol = 1e-10
    end

    @testset "PSDData: variance decomposes across bands" begin
        # Splitting the frequency axis at an interior point (not on a grid node)
        # should give two sub-RMS values whose squares add to psd.rms².
        # Linear interpolation at the split point makes this exact.
        f_mid    = (psd.f_min + psd.f_max) / 2
        rms_low  = rms(psd; f_max = f_mid)
        rms_high = rms(psd; f_min = f_mid)
        @test rms_low < psd.rms
        @test rms_high < psd.rms
        @test rms_low^2 + rms_high^2 ≈ psd.rms^2 rtol = 1e-12
    end

    @testset "PSDData: band outside data range gives zero" begin
        @test rms(psd; f_max = psd.f_min * 0.5) == 0.0
        @test rms(psd; f_min = psd.f_max * 2.0) == 0.0
    end

    dt = 1e-4
    tlist = collect(range(0.0, step = dt, length = 1024))
    trace = AccelerationTrace(
        psd;
        tlist         = tlist,
        f_band_min    = 20.0,
        f_band_max    = 1000.0,
        n_oscillators = 50,
        rng           = Ran2State(Int32(-76989321)),
        target_rms    = psd.rms,
    )
    spectrum = AccelerationSpectrum(trace)

    @testset "AccelerationTrace: returns rms field" begin
        @test rms(trace) === trace.rms
    end

    @testset "AccelerationSpectrum: no windowing matches spectrum.rms" begin
        # The spectrum is normalized so that sum(|F|²)·df = rms²; rms() with
        # default f_min/f_max should recover the same value.
        # The two quantities differ by ~1e-8 in practice: spectrum.rms is the
        # time-domain RMS (sqrt(sum(a²)/N)), whereas rms(spectrum) recomputes
        # the RMS from the FFT output after normalization.  Parseval's theorem
        # guarantees equality in exact arithmetic, but FFT rounding and the
        # normalization factor sqrt(1/df)/N introduce ~1e-8 relative error.
        # rtol=1e-6 is a safe margin above that observed gap.
        @test rms(spectrum) ≈ spectrum.rms rtol = 1e-6
    end

    @testset "AccelerationSpectrum: variance decomposes across bands" begin
        # Splitting the two-sided spectrum at a mid-frequency and summing the
        # squared RMS of each half must reconstruct the total squared RMS.
        # Both ±f sides of each bin must be included via abs(f).
        # We compare to rms(spectrum)² rather than spectrum.rms² because
        # rms_low² + rms_high² is a rearrangement of the exact same discrete
        # sum as rms(spectrum)² = sum(|F|²·df), so they agree to machine
        # precision (rtol=1e-12).  spectrum.rms is the time-domain RMS and
        # differs from rms(spectrum) by the ~1e-8 FFT rounding gap above.
        f_mid = 500.0
        rms_low = rms(spectrum; f_max = f_mid)
        rms_high = rms(spectrum; f_min = f_mid)
        full_rms = rms(spectrum)
        @test rms_low^2 + rms_high^2 ≈ full_rms^2 rtol = 1e-12
    end

    @testset "AccelerationSpectrum: band outside spectrum range gives zero" begin
        # The Nyquist frequency is 0.5/dt = 5000 Hz; asking above that
        # matches no bin and must return exactly 0.
        f_nyquist = 0.5 / dt
        @test rms(spectrum; f_min = f_nyquist + 1.0, f_max = Inf) == 0.0
    end

end


@testset "AccelerationSpectrum reordering" begin

    # Build a pure cosine at a frequency exactly on the FFT grid so there is
    # zero spectral leakage.  With n_samples = 2^10 and dt = 1e-4:
    #   df = f_max / n_samples = 10000 / 1024 Hz
    #   f0 = 10 * df  (10th positive bin)
    # The spectrum should have two and only two non-negligible bins, at ±f0.
    dt        = 1e-4
    n_samples = 2^10
    f_max     = 1.0 / dt
    df        = f_max / n_samples
    f0        = 10 * df   # exactly on-grid: ≈ 97.66 Hz
    A         = 2.0       # amplitude (m/s²)

    tlist        = collect(range(0.0, step = dt, length = n_samples))
    acceleration = A .* cos.(2π * f0 .* tlist)
    displacement = -(A / (2π * f0)^2) .* cos.(2π * f0 .* tlist)
    rms_val      = sqrt(sum(abs2, acceleration) / n_samples)

    trace = AccelerationTrace(
        tlist,
        acceleration,
        displacement,
        zeros(Float64, 3, 1),  # dummy oscillators
        f0,
        f0,                 # dummy band limits
        rms_val,
        n_samples,
        1.0
    )
    spectrum = AccelerationSpectrum(trace)

    @testset "peak frequencies at ±f0" begin
        # For a pure cosine the two dominant spectral bins must be at ±f0.
        power = abs2.(spectrum.spectrum)
        top2 = sortperm(power; rev = true)[1:2]
        top2_f = sort(spectrum.freq_axis[top2])
        @test top2_f[1] ≈ -f0 atol = df / 2
        @test top2_f[2] ≈ f0 atol = df / 2
    end

    @testset "Parseval for pure cosine" begin
        # rms of a cosine A cos(...) is A/√2; the spectrum normalization should
        # reproduce this.
        @test spectrum.rms ≈ A / sqrt(2) rtol = 1e-10
    end

end


@testset "check_acceleration_trace" begin

    psd = PSDData(PSD_FILE; acceleration_unit = 1.0)

    # Generate traces with and without normalization using the same seed.
    # The oscillators drawn are identical in both cases, so raw_rms = trace.rms/η
    # is the same; check_acceleration_trace should return the same n_sigma.
    tlist = collect(range(0.0, step = 1e-4, length = 1024))
    trace_unnorm = AccelerationTrace(
        psd;
        tlist         = tlist,
        f_band_min    = 20.0,
        f_band_max    = 1000.0,
        n_oscillators = 50,
        rng           = Ran2State(Int32(-76989321)),
    )
    trace_norm = AccelerationTrace(
        psd;
        tlist         = tlist,
        f_band_min    = 20.0,
        f_band_max    = 1000.0,
        n_oscillators = 50,
        rng           = Ran2State(Int32(-76989321)),
        target_rms    = psd.rms,
    )

    n_sigma_unnorm = check_acceleration_trace(psd, trace_unnorm)
    n_sigma_norm   = check_acceleration_trace(psd, trace_norm)

    @testset "within 4σ" begin
        # With 50 oscillators the Monte-Carlo scatter is large (O(1/√50)≈14%),
        # but the σ estimate accounts for this.  A correct implementation should
        # produce n_sigma ≤ 4 with probability ~99.994%.  The fixed seed makes
        # this test deterministic.
        @test n_sigma_unnorm ≤ 4.0
        @test n_sigma_norm ≤ 4.0
    end

    @testset "η invariance" begin
        # target_rms only rescales the stored waveform; dividing by η must
        # recover the same raw_rms regardless.  The two n_sigma values should
        # therefore be identical to machine precision.
        @test n_sigma_unnorm ≈ n_sigma_norm rtol = 1e-12
    end

end
