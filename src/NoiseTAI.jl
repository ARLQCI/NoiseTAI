module NoiseTAI

include("ran2.jl")
export Ran2State, ran2!

include("psd.jl")
export PSDData, rms

include("waveform.jl")
export AccelerationTrace, AccelerationSpectrum
export resample_acceleration_trace

include("analysis.jl")
export check_acceleration_trace, check_acceleration_spectrum

include("io.jl")
export write_acceleration_trace, write_acceleration_spectrum

include("main.jl")
export main

end


if abspath(PROGRAM_FILE) == @__FILE__
    NoiseTAI.main()
end
