using PracticalMechanicalSimulation

const MODEL_PATH = joinpath(@__DIR__, "..", "models", "spatial",
    "large-van.lua")
const DURATION = 0.25

allocation_count(timing) = timing.gcstats.malloc + timing.gcstats.realloc +
    timing.gcstats.poolalloc + timing.gcstats.bigalloc

function measure(loaded, repetitions)
    # Compile and exercise both static initialization and dynamics before
    # measuring. Model loading and Lua expansion are intentionally excluded.
    run_spatial_model(loaded; duration = DURATION, samples = 2)

    measurements = NamedTuple[]
    for _ in 1:repetitions
        GC.gc()
        timing = @timed run_spatial_model(loaded;
            duration = DURATION, samples = 2)
        result = timing.value
        statistics = result.solution.stats
        push!(measurements, (;
            time = timing.time,
            allocation_mib = timing.bytes / 2.0^20,
            allocations = allocation_count(timing),
            gc_percent = iszero(timing.time) ? 0.0 :
                100timing.gctime / timing.time,
            static_iterations = result.static_initialization_iterations,
            accepted = statistics.accepted_steps,
            rejected = statistics.rejected_steps,
            residuals = statistics.residual_evaluations,
            jacobians = statistics.jacobian_evaluations,
            newton = statistics.newton_iterations,
            final_norm = maximum(abs, last(result.states))))
    end
    sort!(measurements; by = measurement -> measurement.time)
    measurements[cld(length(measurements), 2)]
end

repetitions = isempty(ARGS) ? 3 : parse(Int, only(ARGS))
repetitions >= 1 || error("the repetition count must be positive")
loaded = load_spatial_model(MODEL_PATH)
result = measure(loaded, repetitions)

println("Large Van: $(length(loaded.active_variable_indices)) active " *
    "variables and $(length(loaded.active_equation_indices)) active equations")
println("Warmed median of $repetitions run(s); model loading is excluded.")
println("Static equilibrium followed by $DURATION simulated seconds:")
println((;
    run_seconds = round(result.time; digits = 6),
    allocated_mib = round(result.allocation_mib; digits = 1),
    result.allocations,
    gc_percent = round(result.gc_percent; digits = 1),
    result.static_iterations,
    accepted_rejected = "$(result.accepted)/$(result.rejected)",
    result.residuals,
    result.jacobians,
    newton_iterations = result.newton,
    final_maximum = result.final_norm))
