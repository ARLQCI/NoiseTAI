using Random: Random


"""
State for the `ran2!` long-period (> 2×10¹⁸) random-number generator.

```julia
rng = Ran2State(-76989321)
```
initialises the generator from *Numerical Recipes in Fortran*, 2nd ed., §7.1
(© 1986–92 Numerical Recipes Software), translated faithfully with `Int32`
arithmetic. Pass a *negative* seed; advance with `ran2!` or `rand(rng)`.

Implements `Random.AbstractRNG` so it can be passed anywhere a Julia RNG is
expected.

# Keyword Arguments
- `verbose`: if `true`, print `"Initialize RAN2"` on the first draw; default `false`
"""
mutable struct Ran2State <: Random.AbstractRNG
    lcg1          :: Int32   # state of first linear congruential generator
    lcg2          :: Int32   # state of second linear congruential generator
    shuffle_table :: Vector{Int32}
    output_reg    :: Int32   # last value drawn from the shuffle table
    verbose       :: Bool
end

Ran2State(seed::Integer; verbose::Bool = false) =
    Ran2State(Int32(seed), Int32(123456789), zeros(Int32, 32), Int32(0), verbose)


"""
Draw the next uniform random `Float64` in (0, 1) from a `Ran2State`.

```julia
x = ran2!(rng)
```
advances `rng` in-place and returns the next value.
"""
function ran2!(state::Ran2State)::Float64
    IM1  = Int32(2147483563)
    IM2  = Int32(2147483399)
    IA1  = Int32(40014)
    IA2  = Int32(40692)
    IQ1  = Int32(53668)
    IQ2  = Int32(52774)
    IR1  = Int32(12211)
    IR2  = Int32(3791)
    NTAB = Int32(32)
    IMM1 = IM1 - Int32(1)
    NDIV = Int32(1) + IMM1 ÷ NTAB
    AM   = 1.0 / Float64(IM1)
    RNMX = 1.0 - 1.2e-7

    if state.lcg1 <= 0
        state.verbose && println("Initialize RAN2")
        state.lcg1 = max(-state.lcg1, Int32(1))
        state.lcg2 = state.lcg1
        for j = (Int(NTAB)+8):-1:1
            q          = state.lcg1 ÷ IQ1
            state.lcg1 = IA1 * (state.lcg1 - q * IQ1) - q * IR1
            if state.lcg1 < 0
                state.lcg1 += IM1
            end
            if j <= Int(NTAB)
                state.shuffle_table[j] = state.lcg1
            end
        end
        state.output_reg = state.shuffle_table[1]
    end

    q          = state.lcg1 ÷ IQ1
    state.lcg1 = IA1 * (state.lcg1 - q * IQ1) - q * IR1
    if state.lcg1 < 0
        state.lcg1 += IM1
    end

    q          = state.lcg2 ÷ IQ2
    state.lcg2 = IA2 * (state.lcg2 - q * IQ2) - q * IR2
    if state.lcg2 < 0
        state.lcg2 += IM2
    end

    table_idx = 1 + Int(state.output_reg ÷ NDIV)
    state.output_reg = state.shuffle_table[table_idx] - state.lcg2
    state.shuffle_table[table_idx] = state.lcg1
    if state.output_reg < 1
        state.output_reg += IMM1
    end

    return min(AM * Float64(state.output_reg), RNMX)
end

Random.rand(rng::Ran2State, ::Random.SamplerTrivial{Random.CloseOpen01{Float64}}) =
    ran2!(rng)
