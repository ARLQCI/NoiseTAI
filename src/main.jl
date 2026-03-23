"""
Run the interactive noise generation program, reading all inputs from stdin.

```julia
main()
```
queries for sampling parameters, a PSD file, oscillator count, random seed,
target RMS, and output file paths, then writes a time trace and spectrum file.
"""
function main()
    println("# noise sample points in time = 2^m. Enter m<=20!")
    log2_n_samples = read_stdin(Int)
    @show log2_n_samples
    n_samples = 2^log2_n_samples
    @show n_samples

    println("tstep (s)?")
    dt = read_stdin(Float64)
    @show dt
    total_time = dt * n_samples
    @show total_time
    f_max = 1.0 / dt
    @show f_max
    df = f_max / n_samples
    @show df

    psd = PSDData()
    show(stdout, MIME"text/plain"(), psd)
    println()

    println("Noise band window (fpmin, fpmax)[Hz]?")
    println("  (Use something inside the above f_min, f_max)")
    f_band_min, f_band_max = read_stdin(Float64, 2)
    @show f_band_min f_band_max

    println("Target number of oscillators nsmax?")
    println("       Do < ntmax = ", n_samples)
    n_oscillators = read_stdin(Int)
    @show n_oscillators

    if n_oscillators > 1_000_000
        @warn "n_oscillators=$n_oscillators exceeds the recommended maximum of 1_000_000"
    end

    println(" dt (s)             : ", dt)
    println(" n_samples          : ", n_samples)
    println(" total_time (s)     : ", total_time)
    println(" f_max (Hz)         : ", f_max)
    println(" Nyquist BW=0.5fmax : ", 0.5 * f_max)
    println(" df (Hz)            : ", df)
    println(" #Sample points     : ", psd.n_samples)
    println(" f_min(sample, Hz)  : ", psd.f_min)
    println(" f_max(sample, Hz)  : ", psd.f_max)
    println(" Noise band min (Hz): ", f_band_min)
    println(" Noise band max (Hz): ", f_band_max)
    println(" # oscillators      : ", n_oscillators)
    println()

    println("Random seed? (int*4) < 0")
    rng_seed = read_stdin(Int32)
    @show rng_seed

    println("Target RMS acceleration, accrms0, in m/s^2 ?")
    println("   Use the above reported arms_sample [m/s^2],")
    println("   or scale it up or down as you wish.")
    target_rms = read_stdin(Float64)
    @show target_rms

    println("Output file for time traces?")
    time_trace_file = read_stdin(String)
    @show time_trace_file
    println("Output file for spectrum?")
    spectrum_file = read_stdin(String)
    @show spectrum_file

    tlist = collect(range(0.0, step = dt, length = n_samples))
    rng = Ran2State(rng_seed; verbose = true)
    trace = AccelerationTrace(
        psd;
        tlist         = tlist,
        f_band_min    = f_band_min,
        f_band_max    = f_band_max,
        n_oscillators = n_oscillators,
        rng           = rng,
        target_rms    = target_rms,
        show_progress = true,
    )

    println()
    show(stdout, MIME"text/plain"(), trace)
    print("\n\n")
    check_acceleration_trace(psd, trace; verbose = true)
    println()

    spectrum = AccelerationSpectrum(trace)
    show(stdout, MIME"text/plain"(), spectrum)
    print("\n\n")
    check_acceleration_spectrum(trace, spectrum; verbose = true)
    println()

    write_acceleration_trace(trace, time_trace_file)
    write_acceleration_spectrum(spectrum, spectrum_file)

    # -- Optional: resample at a different dt (repeating loop) -----------

    while true

        println("Output limited time sample? (1=Yes, else=no)")
        do_resample = 2
        try
            do_resample = read_stdin(Int)
            @show do_resample
        catch
        end
        do_resample != 1 && break

        println("Output file for time traces?")
        resampled_file = read_stdin(String)
        @show resampled_file

        println("tmax(s), dt(s)?")
        resample_duration, resample_dt = read_stdin(Float64, 2)
        @show resample_duration, resample_dt

        n_resample = round(Int, resample_duration / resample_dt)
        tlist_resample = collect(range(0.0, step = resample_dt, length = n_resample + 1))
        resampled =
            resample_acceleration_trace(trace; tlist = tlist_resample, show_progress = true)
        write_acceleration_trace(resampled, resampled_file)

        println(" acc_rms(resampled)             : ", resampled.rms)
        println(" should be close to trace.rms     : ", trace.rms)

    end

end
