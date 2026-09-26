using PracticalMechanicalSimulation

const CASES = [
    (name = "spatial four-bar",
     path = joinpath(@__DIR__, "..", "models", "spatial",
         "torque-driven-four-bar.toml"),
     duration = 1.0),
    (name = "bristle rolling tire",
     path = joinpath(@__DIR__, "..", "models", "spatial",
         "bristle-rolling-tire.toml"),
     duration = 0.1),
]

allocation_count(timing) = timing.gcstats.malloc + timing.gcstats.realloc +
    timing.gcstats.poolalloc + timing.gcstats.bigalloc

function measure(case, repetitions)
    loaded = load_spatial_model(case.path)

    # Warm the complete path before collecting timing or allocation data.
    run_spatial_model(loaded; duration = case.duration, samples = 2)

    measurements = NamedTuple[]
    for _ in 1:repetitions
        GC.gc()
        timing = @timed run_spatial_model(loaded;
            duration = case.duration, samples = 2)
        result = timing.value
        statistics = result.solution.stats
        push!(measurements, (;
            time = timing.time,
            allocation_mib = timing.bytes / 2.0^20,
            allocations = allocation_count(timing),
            gc_percent = iszero(timing.time) ? 0.0 :
                100timing.gctime / timing.time,
            accepted = statistics.accepted_steps,
            rejected = statistics.rejected_steps,
            residuals = statistics.residual_evaluations,
            jacobians = statistics.jacobian_evaluations,
            newton = statistics.newton_iterations,
            final_norm = maximum(abs, last(result.states))))
    end

    # Report one real run rather than combining counters from different runs.
    sort!(measurements; by = measurement -> measurement.time)
    return measurements[cld(length(measurements), 2)]
end

repetitions = isempty(ARGS) ? 3 : parse(Int, only(ARGS))
repetitions >= 1 || error("the repetition count must be positive")

println("Warmed median of $repetitions run(s); model loading is excluded.")
println("| model | run s | allocated MiB | allocations | GC % | " *
        "accepted/rejected | residuals | Jacobians | Newton iterations | " *
        "final maximum |")
println("|:--|--:|--:|--:|--:|:--|--:|--:|--:|--:|")
for case in CASES
    result = measure(case, repetitions)
    println("| ", case.name,
        " | ", round(result.time; digits = 6),
        " | ", round(result.allocation_mib; digits = 1),
        " | ", result.allocations,
        " | ", round(result.gc_percent; digits = 1),
        " | ", result.accepted, "/", result.rejected,
        " | ", result.residuals,
        " | ", result.jacobians,
        " | ", result.newton,
        " | ", result.final_norm,
        " |")
end
