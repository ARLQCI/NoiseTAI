using DataInterpolations: DataInterpolations, LinearInterpolation, PCHIPInterpolation
using RecipesBase: RecipesBase, @recipe


"""
Acceleration PSD S(f) loaded from a two-column ASCII file.

```julia
psd = PSDData(file; acceleration_unit=1.0, interpolation_method=:linear)
```
reads the power spectral density (PSD) from `file`. The `file` must contain
two whitespace-separated columns; lines
starting with `#` are ignored as comments:

* Column 1: frequency (Hz), must be sorted in ascending order
* Column 2: one-sided acceleration PSD in a²/Hz

Calling `PSDData()` without positional arguments instead
prompts for the file path and `acceleration_unit` interactively from stdin.

The `acceleration_unit` is the conversion factor from the file's acceleration
unit ``a`` to m/s²: `acceleration_unit=9.81` if the PSD is in g²/Hz;
`acceleration_unit=1` if already in (m/s²)²/Hz.

By the Wiener–Khinchin theorem, √∫ S(f) df is the RMS √⟨a²⟩ of a stationary
random process with one-sided PSD S(f), so `psd.rms` is the physically
expected RMS acceleration.

The resulting `psd` object is a callable `psd(f)` → S(f) in (m/s²)²/Hz,
regardless of the file's unit, returning interpolated values for the PSD.

A `PSDData` object can be plotted directly with any `Plots`-compatible backend,
showing S(f) in (m/s²)²/Hz vs frequency in Hz.

# Keyword Arguments
- `acceleration_unit`: conversion factor from file units to m/s²; default `1.0`
- `interpolation_method`: `:linear` (default, reproduces trapezoid-rule results
  exactly) or `:pchip` (monotone-preserving cubic, smoother S(f) between nodes)

# Fields
- `f_min`: lowest frequency in the PSD file (Hz)
- `f_max`: highest frequency in the PSD file (Hz)
- `n_samples`: number of data points in the PSD file
- `rms`: ``√∫S(f)df`` in (m/s²), the expected root-mean-square of the
  acceleration trace described by the  PSD.
"""
struct PSDData{Itp}
    f_min::Float64
    f_max::Float64
    n_samples::Int
    rms::Float64
    itp::Itp  # interpolant (internal)
end


function PSDData(
    file::String;
    acceleration_unit::Float64 = 1.0,
    interpolation_method::Symbol = :linear
)
    freqs  = Float64[]
    values = Float64[]
    open(file, "r") do fh
        for line in eachline(fh)
            s = strip(line)
            (isempty(s) || startswith(s, '#')) && continue
            cols = split(s)
            length(cols) < 2 && continue
            push!(freqs, parse(Float64, replace(cols[1], "d" => "e")))
            push!(values, parse(Float64, replace(cols[2], "d" => "e")))
        end
    end
    values .*= acceleration_unit^2  # convert to (m/s²)²/Hz

    itp = if interpolation_method == :linear
        LinearInterpolation(values, freqs)
    elseif interpolation_method == :pchip
        PCHIPInterpolation(values, freqs)
    else
        error("Unknown interpolation_method: $interpolation_method. Use :linear or :pchip.")
    end

    rms_sq = DataInterpolations.integral(itp, freqs[1], freqs[end])
    return PSDData(freqs[1], freqs[end], length(freqs), sqrt(rms_sq), itp)
end


function PSDData(; interpolation_method::Symbol = :linear)
    gravity = 9.81
    println("Acceleration PSD sample file?")
    println(" (MUST BE sorted in ascending frequency)")
    psd_file = read_stdin(String)
    @show psd_file

    println("Acceleration Unit G (m/s^2) = ", gravity)
    println("What is the unit for a, in m/s^2, in your sample?")
    acceleration_unit = read_stdin(Float64)
    @show acceleration_unit

    return PSDData(
        psd_file;
        acceleration_unit = acceleration_unit,
        interpolation_method = interpolation_method
    )
end


function (psd::PSDData)(f::Real)
    if f < psd.f_min || f > psd.f_max
        error("Frequency $f Hz is outside PSD data range [$(psd.f_min), $(psd.f_max)] Hz")
    end
    return psd.itp(f)
end


function Base.show(io::IO, psd::PSDData)
    print(io, "PSDData(f_min=", psd.f_min, ", f_max=", psd.f_max, ", rms=", psd.rms, ")")
end


function Base.show(io::IO, ::MIME"text/plain", psd::PSDData)
    println(io, "PSDData S(f) with ", psd.n_samples, " samples")
    println(io, "   f_min: ", psd.f_min, " Hz")
    println(io, "   f_max: ", psd.f_max, " Hz")
    print(io, "   rms:  ", psd.rms, " m/s²")
end


@recipe function f(psd::PSDData)
    xguide --> "Frequency (Hz)"
    yguide --> "S(f) ((m/s²)²/Hz)"
    psd.itp.t, psd.itp.u
end


"""
RMS acceleration from a `PSDData`, `AccelerationTrace`, or `AccelerationSpectrum`.

```julia
band_rms  = rms(psd; f_min=20.0, f_max=1000.0)
trace_rms = rms(trace)
spec_rms  = rms(spectrum; f_min=20.0, f_max=1000.0)
```

For [`PSDData`](@ref), integrates S(f) over [f_min, f_max] using the interpolant;
without windowing (`f_min=0`, `f_max=Inf`), returns `psd.rms`. For
[`AccelerationTrace`](@ref), returns `trace.rms` directly (no keyword arguments).
For [`AccelerationSpectrum`](@ref), sums |F(f)|²·df over all f with
|f| ∈ [f_min, f_max]; without windowing, recovers `spectrum.rms`.

Comparing `rms(psd; f_min=f_band_min, f_max=f_band_max)` to `rms(trace)` when
`target_rms=nothing` isolates Monte-Carlo fluctuations from the systematic
band-limitation bias.

# Keyword Arguments
- `f_min`: lower frequency bound (Hz); default `0.0`
- `f_max`: upper frequency bound (Hz); default `Inf`
"""
function rms(psd::PSDData; f_min::Float64 = 0.0, f_max::Float64 = Inf)
    f_min == 0.0 && f_max == Inf && return psd.rms
    a = max(f_min, psd.f_min)
    b = min(f_max, psd.f_max)
    b <= a && return 0.0
    return sqrt(DataInterpolations.integral(psd.itp, a, b))
end
