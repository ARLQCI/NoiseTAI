using Printf: Printf, @printf


"""
Read one stdin line and parse `n` values of the given type.

```julia
(f_min, f_max) = read_stdin(Float64, 2)
x              = read_stdin(Int)
```
strips any label up to and including the first colon, so `"tstep (s): 1e-4"`
and `"1e-4"` parse identically. `Float64` inputs accept Fortran syntax
(`"1d-4"` → `1e-4`). Returns a `Tuple` for `n > 1`, a scalar for `n = 1`.
"""
function read_stdin(type, n)
    line = strip(readline())
    text = strip(replace(line, r"^[^:]*:" => ""))
    data = split(text)
    values = similar(data, type)
    for (i, item) in enumerate(data)
        if type == Float64
            item = replace(item, "d" => "e")
        end
        values[i] = type == String ? String(item) : parse(type, item)
    end
    return Tuple(values)
end

read_stdin(type) = read_stdin(type, 1)[1]  # single value


"""
Write an [`AccelerationTrace`](@ref) to a time-trace ASCII file.

```julia
write_acceleration_trace(trace, "time_trace.dat")
```
writes three columns: time (s), acceleration (m/s²), displacement (m).
The time column is taken directly from `trace.tlist`.
"""
function write_acceleration_trace(trace::AccelerationTrace, time_trace_file::String)
    open(time_trace_file, "w") do fh
        println(fh, "# Acceleration time trace")
        println(fh, "# Column 1: time (s)")
        println(fh, "# Column 2: acceleration (m/s^2)")
        println(fh, "# Column 3: displacement (m)")
        for (t, a, d) in zip(trace.tlist, trace.acceleration, trace.displacement)
            @printf(fh, "%15.8E %15.8E %15.8E\n", t, a, d)
        end
    end
end


"""
Write an [`AccelerationSpectrum`](@ref) to a spectrum ASCII file.

```julia
write_acceleration_spectrum(spectrum, "spectrum.dat")
```
writes four columns: frequency (Hz), Re (m/s²/√Hz), Im (m/s²/√Hz),
two-sided PSD ((m/s²)²/Hz).
"""
function write_acceleration_spectrum(spectrum::AccelerationSpectrum, spectrum_file::String)
    open(spectrum_file, "w") do fh
        println(fh, "# Two-sided acceleration spectrum")
        println(fh, "# Normalized so that integral of PSD over all f equals RMS^2")
        println(fh, "# Column 1: frequency (Hz)")
        println(fh, "# Column 2: real part of spectral amplitude (m/s^2/sqrt(Hz))")
        println(fh, "# Column 3: imaginary part of spectral amplitude (m/s^2/sqrt(Hz))")
        println(fh, "# Column 4: two-sided power spectral density (m/s^2)^2/Hz")
        for (f, s) in zip(spectrum.freq_axis, spectrum.spectrum)
            @printf(fh, "%15.8E %15.8E %15.8E %15.8E\n", f, real(s), imag(s), abs2(s))
        end
    end
end
