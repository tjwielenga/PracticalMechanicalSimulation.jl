#!/usr/bin/env julia

using Dates
using LinearAlgebra
using Printf
using SHA
using Statistics
using TOML
using PracticalMechanicalSimulation

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RELEASE = "v0.2.0"
const RELEASE_COMMIT = "b7189d3eb62ca99b27facae78e7d0e6202ac2a62"
const OUTPUT_DIRECTORY = joinpath(@__DIR__, RELEASE[2:end])

include(joinpath(ROOT, "benchmark", "PendulumChainBenchmark.jl"))
using .PendulumChainBenchmark

module BladeEvidence
include(joinpath(@__DIR__, "..", "..", "examples", "spatial",
    "rotating_flexible_blade.jl"))
end

function parse_arguments(arguments)
    repetitions = 3
    quick = false
    for argument in arguments
        if argument == "--quick"
            quick = true
        elseif startswith(argument, "--repeats=")
            repetitions = parse(Int, split(argument, "="; limit = 2)[2])
        else
            error("unknown argument: $argument")
        end
    end
    repetitions >= 1 || error("the repetition count must be positive")
    quick && (repetitions = 1)
    (; repetitions, quick)
end

function git_output(arguments...)
    strip(read(`git -C $ROOT $(arguments)`, String))
end

function release_source_is_unchanged()
    process = run(ignorestatus(`git -C $ROOT diff --quiet $RELEASE_COMMIT -- src Project.toml`))
    success(process)
end

function allocation_count(timing)
    timing.gcstats.malloc + timing.gcstats.realloc +
        timing.gcstats.poolalloc + timing.gcstats.bigalloc
end

function median_time_run(measure, repetitions)
    rows = [measure() for _ in 1:repetitions]
    sort!(rows; by = row -> row.run_seconds)
    rows[cld(length(rows), 2)]
end

function csv_field(value)
    text = value isa AbstractFloat ? repr(value) : string(value)
    if occursin(r"[\",\n\r]", text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    text
end

function write_csv(path, rows)
    isempty(rows) && error("cannot write an empty CSV table")
    columns = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((csv_field(getproperty(row, column))
                for column in columns), ','))
        end
    end
end

function write_metadata(options)
    cpu = first(Sys.cpu_info())
    project = TOML.parsefile(joinpath(ROOT, "Project.toml"))
    paper_manifest = joinpath(ROOT, "test", "paper", "Manifest.toml")
    dirty = split(git_output("status", "--short"), '\n'; keepempty = false)
    metadata = Dict(
        "evidence" => Dict(
            "generated_utc" => string(now(UTC)),
            "release" => RELEASE,
            "release_commit" => RELEASE_COMMIT,
            "checkout_commit" => git_output("rev-parse", "HEAD"),
            "simulator_source_matches_release" =>
                release_source_is_unchanged(),
            "working_tree_changes" => dirty,
            "timed_repetitions" => options.repetitions,
            "quick_run" => options.quick,
            "paper_environment_manifest_sha256" => isfile(paper_manifest) ?
                bytes2hex(sha256(read(paper_manifest))) : "missing",
        ),
        "software" => Dict(
            "package_version" => project["version"],
            "julia_version" => string(VERSION),
            "blas" => string(BLAS.get_config()),
            "julia_threads" => Threads.nthreads(),
        ),
        "computer" => Dict(
            "kernel" => string(Sys.KERNEL),
            "architecture" => string(Sys.ARCH),
            "cpu_model" => cpu.model,
            "logical_cores" => length(Sys.cpu_info()),
            "word_size" => Sys.WORD_SIZE,
        ),
        "protocol" => Dict(
            "timing" => "median-time complete run after path warm-up",
            "open_and_closed_chain_duration_s" => 0.25,
            "rotor_train_duration_s" => 1.0,
            "bouncing_ball_duration_s" => 1.5,
            "large_van_microbenchmark_dynamic_duration_s" => 0.25,
            "large_van_maneuver_dynamic_duration_s" => 10.0,
            "relative_tolerance" => 1.0e-7,
            "absolute_tolerance" => 1.0e-9,
        ),
    )
    open(joinpath(OUTPUT_DIRECTORY, "metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted = true)
    end
    metadata
end

function compact_benchmark_row(row, size_name)
    (; size = getproperty(row, size_name),
       variables = row.variables,
       states = row.states,
       jacobian_nonzeros = row.jacobian_nonzeros,
       load_seconds = row.load_seconds,
       run_seconds = row.run_seconds,
       allocated_mib = row.run_mebibytes,
       accepted_steps = row.accepted_steps,
       rejected_steps = row.rejected_steps,
       corrector_failures = row.corrector_failures,
       residual_evaluations = row.residual_evaluations,
       jacobian_evaluations = row.jacobian_evaluations,
       numerical_factorizations = row.factorizations,
       symbolic_factorizations = row.symbolic_factorizations,
       newton_iterations = row.newton_iterations,
       maximum_order = row.maximum_order,
       maximum_joint_error = row.maximum_joint_error)
end

function measure_scalable_case(case, warm, counts, repetitions, size_name)
    warm()
    rows = NamedTuple[]
    for count in counts
        println("Measuring $size_name = $count")
        row = median_time_run(() -> case(count), repetitions)
        push!(rows, compact_benchmark_row(row, size_name))
    end
    rows
end

function measure_open_chains(counts, repetitions)
    measure_scalable_case(
        count -> PendulumChainBenchmark.benchmark_case(count;
            duration = 0.25, samples = 11),
        () -> run_pendulum_chain(2; duration = 0.02, samples = 2),
        counts, repetitions, :links)
end

function measure_closed_chains(counts, repetitions)
    measure_scalable_case(
        count -> PendulumChainBenchmark.parallelogram_benchmark_case(count;
            duration = 0.25, samples = 11),
        () -> run_parallelogram_chain(1; duration = 0.02, samples = 2),
        counts, repetitions, :cells)
end

function measure_rotor_trains(counts, repetitions)
    run_rotor_train(2; duration = 0.02, samples = 3)
    rows = NamedTuple[]
    for count in counts
        println("Measuring rotors = $count")
        result = median_time_run(() ->
            PendulumChainBenchmark.rotor_benchmark_case(count;
                duration = 1.0, samples = 201), repetitions)
        push!(rows, (;
            rotors = count,
            variables = result.variables,
            states = result.states,
            run_seconds = result.run_seconds,
            allocated_mib = result.run_mebibytes,
            accepted_steps = result.accepted_steps,
            rejected_steps = result.rejected_steps,
            corrector_failures = result.corrector_failures,
            numerical_factorizations = result.factorizations,
            symbolic_factorizations = result.symbolic_factorizations,
            minimum_frequency_rad_s = result.minimum_frequency,
            maximum_frequency_rad_s = result.maximum_frequency,
            maximum_angle_error_rad = result.maximum_angle_error,
            maximum_angular_velocity_error_rad_s =
                result.maximum_angular_velocity_error,
            final_energy_ratio = result.final_energy_ratio))
    end
    rows
end

function damped_chain_source(choice)
    count = 10
    digits = 2
    link_name(index) = "link" * lpad(string(index), digits, '0')
    joint_name(index) = "joint" * lpad(string(index), digits, '0')
    source = pendulum_chain_toml(count;
        initial_angle = pi / 2, end_time = 10.0, output_samples = 1001)
    io = IOBuffer()
    println(io, source)
    for index in 1:count
        first_marker = "$(link_name(index)).top"
        second_marker = index == 1 ? "ground.origin" :
            "$(link_name(index - 1)).bottom"
        println(io, "\n[$(joint_name(index))_damper]")
        println(io, "type = \"torsional_spring_damper\"")
        println(io, "markers = [\"$first_marker\", \"$second_marker\"]")
        println(io, "stiffness = 0.0")
        println(io, "damping = 0.01")
    end
    if choice == :body
        println(io, "\n[state_selection]")
        println(io, "method = \"preferred\"")
        values = ["\"$(link_name(index)).omega\"" for index in 1:count]
        println(io, "preferred_velocities = [", join(values, ", "), "]")
        println(io, "allow_fallback = false")
    elseif choice == :relative
        source_with_dampers = String(take!(io))
        for index in 1:count
            joint = joint_name(index)
            pattern = Regex("(\\[$joint\\]\\ntype = \\\"revolute\\\"\\n" *
                "markers = \\[[^\\n]+\\]\\n)")
            source_with_dampers = replace(source_with_dampers,
                pattern => s"\1rotation_coordinates = true\n")
        end
        io = IOBuffer()
        println(io, source_with_dampers)
        println(io, "\n[state_selection]")
        println(io, "method = \"preferred\"")
        values = ["\"$(joint_name(index)).omega\"" for index in 1:count]
        println(io, "preferred_velocities = [", join(values, ", "), "]")
        println(io, "allow_fallback = false")
    elseif choice != :automatic
        error("unknown state choice: $choice")
    end
    String(take!(io))
end

function measure_state_choices(repetitions)
    choices = ((:automatic, "Automatic QR"),
        (:body, "Preferred body angular velocities"),
        (:relative, "Preferred relative joint angular velocities"))
    sources = Dict(choice => damped_chain_source(choice)
        for (choice, _) in choices)
    for (choice, _) in choices
        run_planar_model(IOBuffer(sources[choice]))
    end
    rows = NamedTuple[]
    for (choice, label) in choices
        println("Measuring state choice: $label")
        measurement = median_time_run(function ()
            GC.gc()
            timing = @timed run_planar_model(IOBuffer(sources[choice]))
            result = timing.value
            stats = result.solution.stats
            (; run_seconds = timing.time,
               allocated_mib = timing.bytes / 2.0^20,
               variables = length(result.loaded.layout.catalog.variables),
               selected_states = result.loaded.analysis.degrees_of_freedom,
               accepted_steps = stats.accepted_steps,
               rejected_steps = stats.rejected_steps,
               residual_evaluations = stats.residual_evaluations,
               numerical_factorizations = stats.factorizations,
               symbolic_factorizations = stats.symbolic_factorizations,
               newton_iterations = stats.newton_iterations,
               corrector_failures = stats.corrector_failures)
        end, repetitions)
        push!(rows, merge((; state_choice = label), measurement))
    end
    rows
end

function measure_bouncing_balls(counts, repetitions)
    run_bouncing_ball_bank(1; duration = 1.5, samples = 5)
    rows = NamedTuple[]
    for count in counts
        println("Measuring soft contact restart, balls = $count")
        result = median_time_run(() ->
            PendulumChainBenchmark.contact_benchmark_case(count;
                duration = 1.5, samples = 301), repetitions)
        push!(rows, (;
            balls = count,
            policy = "soft",
            variables = result.variables,
            run_seconds = result.run_seconds,
            allocated_mib = result.run_mebibytes,
            accepted_steps = result.accepted_steps,
            rejected_steps = result.rejected_steps,
            corrector_failures = result.corrector_failures,
            root_evaluations = result.root_evaluations,
            events = result.events_found,
            history_restarts = result.history_restarts,
            maximum_sampled_penetration_m =
                result.maximum_sampled_penetration))
    end
    rows
end

function measure_modal_verification()
    cases = (
        ("Planar torsional pendulum", :planar,
            joinpath(ROOT, "models", "planar", "modal-pendulum.toml")),
        ("Spatial torsional pendulum", :spatial,
            joinpath(ROOT, "models", "spatial",
                "modal-revolute-pendulum.toml")),
        ("Spatial three-link pendulum", :spatial,
            joinpath(ROOT, "models", "spatial",
                "modal-three-link-pendulum.toml")),
    )
    rows = NamedTuple[]
    for (name, dimension, path) in cases
        println("Measuring modal verification: $name")
        result = dimension == :planar ? run_planar_model(path) :
            run_spatial_model(path)
        reference = occursin("three-link", name) ?
            "independent state matrix" : "closed form"
        for mode in eachindex(result.natural_frequencies_hz)
            push!(rows, (;
                model = name,
                mode,
                natural_frequency_hz = result.natural_frequencies_hz[mode],
                damped_frequency_hz = result.damped_frequencies_hz[mode],
                damping_ratio = result.damping_ratios[mode],
                reference,
                equation_error = result.equation_errors[mode]))
        end
    end
    rows
end

function measure_blades(quick)
    cases = quick ? ((4, 7.5),) : ((4, 0.0), (4, 7.5), (8, 7.5))
    rows = NamedTuple[]
    for (segments, speed) in cases
        println("Measuring rotating blade: $segments segments at $speed rad/s")
        result = BladeEvidence.run_case(segments, speed)
        stats = result.dynamic.solution.stats
        push!(rows, (;
            segments,
            target_speed_rad_s = speed,
            actual_speed_rad_s = result.speed,
            static_tip_height_m = result.initial_tip_z,
            final_tip_height_m = result.final_tip_z,
            first_frequency_hz = result.first_frequency,
            maximum_modal_equation_error = result.equation_error,
            accepted_steps = stats.accepted_steps,
            rejected_steps = stats.rejected_steps,
            corrector_failures = stats.corrector_failures,
            state_reselections = stats.state_reselections))
    end
    rows
end

function measure_van(repetitions)
    path = joinpath(ROOT, "models", "spatial", "large-van.lua")
    loaded = load_spatial_model(path)
    run_spatial_model(loaded; duration = 0.25, samples = 2)
    median_time_run(function ()
        GC.gc()
        timing = @timed run_spatial_model(loaded;
            duration = 0.25, samples = 2)
        result = timing.value
        stats = result.solution.stats
        (; implementation = "0.2.0 current implementation",
           active_variables = length(loaded.active_variable_indices),
           run_seconds = timing.time,
           allocated_mib = timing.bytes / 2.0^20,
           allocations = allocation_count(timing),
           gc_seconds = timing.gctime,
           static_iterations = result.static_initialization_iterations,
           accepted_steps = stats.accepted_steps,
           rejected_steps = stats.rejected_steps,
           residual_evaluations = stats.residual_evaluations,
           jacobian_evaluations = stats.jacobian_evaluations,
           numerical_factorizations = stats.factorizations,
           symbolic_factorizations = stats.symbolic_factorizations,
           newton_iterations = stats.newton_iterations,
           final_maximum = maximum(abs, last(result.states)))
    end, repetitions)
end

function measure_van_maneuver(repetitions)
    path = joinpath(ROOT, "models", "spatial", "large-van.lua")
    loaded = load_spatial_model(path)
    body = loaded.bodies[Symbol("van.body")]
    rotation_matrix =
        PracticalMechanicalSimulation.SpatialComponentAssembly.rotation_matrix

    # Compile the same static-initialization and dynamic path before timing.
    run_spatial_model(loaded; duration = 0.25, samples = 2)
    median_time_run(function ()
        GC.gc()
        timing = @timed run_spatial_model(loaded;
            duration = 10.0, samples = 601)
        result = timing.value
        stats = result.solution.stats
        ground_speeds = [norm(state[body.velocity_variables])
            for state in result.states]
        minimum_index = argmin(ground_speeds)
        initial_state = first(result.states)
        final_state = last(result.states)
        initial_body_velocity = transpose(rotation_matrix(
            initial_state[body.euler_parameter_variables])) *
            initial_state[body.velocity_variables]
        final_body_velocity = transpose(rotation_matrix(
            final_state[body.euler_parameter_variables])) *
            final_state[body.velocity_variables]
        (; maneuver = "180 degree steering-wheel sweep",
           simulated_seconds = 10.0,
           output_samples = length(result.times),
           run_seconds = timing.time,
           simulated_seconds_per_run_second = 10.0 / timing.time,
           allocated_mib = timing.bytes / 2.0^20,
           gc_seconds = timing.gctime,
           static_iterations = result.static_initialization_iterations,
           initial_ground_speed_m_s = first(ground_speeds),
           minimum_ground_speed_m_s = ground_speeds[minimum_index],
           minimum_speed_time_s = result.times[minimum_index],
           final_ground_speed_m_s = last(ground_speeds),
           initial_body_forward_speed_m_s = -initial_body_velocity[1],
           final_body_forward_speed_m_s = -final_body_velocity[1],
           accepted_steps = stats.accepted_steps,
           rejected_steps = stats.rejected_steps,
           residual_evaluations = stats.residual_evaluations,
           jacobian_evaluations = stats.jacobian_evaluations,
           numerical_factorizations = stats.factorizations,
           symbolic_factorizations = stats.symbolic_factorizations,
           newton_iterations = stats.newton_iterations,
           corrector_failures = stats.corrector_failures,
           state_reselections = stats.state_reselections)
    end, repetitions)
end

function markdown_table(io, rows, columns)
    headers = first.(columns)
    fields = last.(columns)
    println(io, "| ", join(headers, " | "), " |")
    println(io, "|", join(fill(":--", length(columns)), "|"), "|")
    for row in rows
        values = [markdown_value(getproperty(row, field)) for field in fields]
        println(io, "| ", join(values, " | "), " |")
    end
end

function markdown_value(value)
    value isa AbstractFloat || return string(value)
    iszero(value) && return "0"
    magnitude = abs(value)
    if magnitude < 1.0e-4 || magnitude >= 1.0e5
        return @sprintf("%.4e", value)
    end
    string(round(value; sigdigits = 6))
end

function write_tables(all_rows)
    open(joinpath(OUTPUT_DIRECTORY, "tables.md"), "w") do io
        println(io, "# Tables generated for $RELEASE\n")
        println(io, "Unrounded source values are in the adjacent CSV files.\n")
        for (title, rows) in all_rows
            println(io, "## ", title, "\n")
            columns = [(replace(string(key), '_' => ' '), key)
                for key in keys(first(rows))]
            markdown_table(io, rows, columns)
            println(io)
        end
    end
end

function write_scope()
    open(joinpath(OUTPUT_DIRECTORY, "scope.md"), "w") do io
        println(io, "# Evidence scope for $RELEASE\n")
        println(io, "The CSV files in this directory are new measurements of " *
            "the released 0.2.0 simulation program.\n")
        println(io, "The following comparisons are intentionally not regenerated:")
        println(io, "\n- Table 20's generic sparse LU, first UMFPACK, and " *
            "symbolic-only rows describe earlier implementations. The open-chain " *
            "CSV supplies the current final implementation.")
        println(io, "- Table 25's hard-restart rows describe the former restart " *
            "policy. Version 0.2.0 uses the soft restart; the CSV measures it.")
        println(io, "- Table 26's ordinary-array and first fixed-size rows describe " *
            "earlier implementations. `large_van.csv` supplies the current final row.")
        println(io, "\nTables 3--5 are regenerated by " *
            "`formulation_evidence.jl` in the isolated `test/paper` environment " *
            "because they compare historical formulations rather than the " *
            "supported model runner. Their values are stored in the three " *
            "`formulation_*.csv` files.")
    end
end

function main(arguments = ARGS)
    options = parse_arguments(arguments)
    mkpath(OUTPUT_DIRECTORY)
    metadata = write_metadata(options)
    metadata["evidence"]["simulator_source_matches_release"] || error(
        "src/ or Project.toml differs from release commit $RELEASE_COMMIT")

    paper_project = joinpath(ROOT, "test", "paper")
    formulation_exporter = joinpath(@__DIR__, "formulation_evidence.jl")
    println("Measuring isolated historical formulations")
    run(`$(Base.julia_cmd()) --project=$paper_project $formulation_exporter`)

    counts = options.quick ? (10,) : (10, 25, 50)
    tables = Pair{String,Any}[]

    modal = measure_modal_verification()
    write_csv(joinpath(OUTPUT_DIRECTORY, "modal_verification.csv"), modal)
    push!(tables, "Modal verification" => modal)

    open_chains = measure_open_chains(counts, options.repetitions)
    write_csv(joinpath(OUTPUT_DIRECTORY, "open_pendulum_chains.csv"), open_chains)
    push!(tables, "Open pendulum chains" => open_chains)

    state_choices = measure_state_choices(options.repetitions)
    write_csv(joinpath(OUTPUT_DIRECTORY, "state_choices.csv"), state_choices)
    push!(tables, "State choices" => state_choices)

    closed_chains = measure_closed_chains(counts, options.repetitions)
    write_csv(joinpath(OUTPUT_DIRECTORY, "closed_parallelogram_chains.csv"),
        closed_chains)
    push!(tables, "Closed parallelogram chains" => closed_chains)

    rotors = measure_rotor_trains(counts, options.repetitions)
    write_csv(joinpath(OUTPUT_DIRECTORY, "rotor_trains.csv"), rotors)
    push!(tables, "Rotor trains" => rotors)

    blades = measure_blades(options.quick)
    write_csv(joinpath(OUTPUT_DIRECTORY, "rotating_flexible_blades.csv"), blades)
    push!(tables, "Rotating flexible blades" => blades)

    contacts = measure_bouncing_balls(counts, options.repetitions)
    write_csv(joinpath(OUTPUT_DIRECTORY, "bouncing_ball_soft_restarts.csv"),
        contacts)
    push!(tables, "Bouncing-ball soft restarts" => contacts)

    van = [measure_van(options.repetitions)]
    write_csv(joinpath(OUTPUT_DIRECTORY, "large_van.csv"), van)
    push!(tables, "Large Van implementation microbenchmark" => van)

    van_maneuver = [measure_van_maneuver(options.repetitions)]
    write_csv(joinpath(OUTPUT_DIRECTORY, "large_van_10s_maneuver.csv"),
        van_maneuver)
    push!(tables, "Large Van 10-second maneuver" => van_maneuver)

    write_tables(tables)
    write_scope()
    println("Evidence written to $OUTPUT_DIRECTORY")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
