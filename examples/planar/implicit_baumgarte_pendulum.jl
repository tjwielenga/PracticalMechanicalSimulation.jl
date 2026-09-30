using LinearAlgebra
using OrdinaryDiffEq
import DASSL
import Sundials

if !isdefined(@__MODULE__, :HistoricalDDASSL)
    include(joinpath(@__DIR__, "..", "..", "src", "common", "HistoricalDDASSL.jl"))
end

if !isdefined(@__MODULE__, :implicit_initial_conditions)
    include(joinpath(@__DIR__, "implicit_displacement_pendulum.jl"))
end

struct BaumgarteParameters{P,T}
    mechanical::P
    correction_time::T
end

correction_frequency(p::BaumgarteParameters) = inv(p.correction_time)

"""Critically damped acceleration-level constraint equation."""
function implicit_baumgarte_pendulum!(equations, dy, y,
                                       p::BaumgarteParameters, t)
    mechanical = p.mechanical
    frequency = correction_frequency(p)
    r_global, d_global = marker_vectors(y[3], mechanical)

    equations[1:2] .= mechanical.mass .* dy[4:5] .- y[7:8] .-
                      mechanical.mass .* mechanical.gravity
    equations[3] = mechanical.inertia * dy[6] - dot(d_global, y[7:8])
    equations[4:5] .= dy[1:2] .- y[4:5]
    equations[6] = dy[3] - y[6]

    phi = y[1:2] + r_global - mechanical.pin
    phi_dot = y[4:5] + d_global .* y[6]
    phi_ddot = dy[4:5] + d_global .* dy[6] - r_global .* y[6]^2
    equations[7:8] .= phi_ddot .+ 2frequency .* phi_dot .+
                       frequency^2 .* phi
    return nothing
end

"""Analytical matrix F_y + coefficient*F_dy."""
function implicit_baumgarte_jacobian(t, y, dy, coefficient,
                                      p::BaumgarteParameters)
    mechanical = p.mechanical
    frequency = correction_frequency(p)
    r_global, d_global = marker_vectors(y[3], mechanical)
    jacobian = zeros(eltype(y), 8, 8)

    jacobian[1, 4] = coefficient * mechanical.mass
    jacobian[2, 5] = coefficient * mechanical.mass
    jacobian[1, 7] = -1
    jacobian[2, 8] = -1
    jacobian[3, 3] = dot(r_global, y[7:8])
    jacobian[3, 6] = coefficient * mechanical.inertia
    jacobian[3, 7:8] .= -d_global

    jacobian[4, 1] = coefficient
    jacobian[5, 2] = coefficient
    jacobian[4, 4] = -1
    jacobian[5, 5] = -1
    jacobian[6, 3] = coefficient
    jacobian[6, 6] = -1

    jacobian[7, 1] = frequency^2
    jacobian[8, 2] = frequency^2
    jacobian[7:8, 3] .= -r_global .* dy[6] .-
        d_global .* y[6]^2 .- 2frequency .* r_global .* y[6] .+
        frequency^2 .* d_global
    jacobian[7, 4] = coefficient + 2frequency
    jacobian[8, 5] = coefficient + 2frequency
    jacobian[7:8, 6] .= coefficient .* d_global .-
        2 .* r_global .* y[6] .+ 2frequency .* d_global
    return jacobian
end

function implicit_baumgarte_jacobian!(jacobian, dy, y,
                                       p::BaumgarteParameters,
                                       coefficient, t)
    jacobian .= implicit_baumgarte_jacobian(t, y, dy, coefficient, p)
    return nothing
end

function historical_baumgarte_residual!(equations, t, y, dy,
                                         p::BaumgarteParameters)
    implicit_baumgarte_pendulum!(equations, dy, y, p, t)
    return nothing
end

function historical_baumgarte_jacobian!(jacobian, t, y, dy, coefficient,
                                         p::BaumgarteParameters)
    jacobian .= implicit_baumgarte_jacobian(t, y, dy, coefficient, p)
    return nothing
end

"""First-order stabilization using the velocity and displacement constraints."""
function implicit_first_order_stabilized_pendulum!(equations, dy, y,
                                                    p::BaumgarteParameters, t)
    mechanical = p.mechanical
    frequency = correction_frequency(p)
    r_global, d_global = marker_vectors(y[3], mechanical)

    equations[1:2] .= mechanical.mass .* dy[4:5] .- y[7:8] .-
                      mechanical.mass .* mechanical.gravity
    equations[3] = mechanical.inertia * dy[6] - dot(d_global, y[7:8])
    equations[4:5] .= dy[1:2] .- y[4:5]
    equations[6] = dy[3] - y[6]

    phi = y[1:2] + r_global - mechanical.pin
    phi_dot = y[4:5] + d_global .* y[6]
    equations[7:8] .= phi_dot .+ frequency .* phi
    return nothing
end

"""Analytical matrix for first-order stabilization."""
function implicit_first_order_stabilized_jacobian(t, y, dy, coefficient,
                                                   p::BaumgarteParameters)
    mechanical = p.mechanical
    frequency = correction_frequency(p)
    r_global, d_global = marker_vectors(y[3], mechanical)
    jacobian = zeros(eltype(y), 8, 8)

    jacobian[1, 4] = coefficient * mechanical.mass
    jacobian[2, 5] = coefficient * mechanical.mass
    jacobian[1, 7] = -1
    jacobian[2, 8] = -1
    jacobian[3, 3] = dot(r_global, y[7:8])
    jacobian[3, 6] = coefficient * mechanical.inertia
    jacobian[3, 7:8] .= -d_global

    jacobian[4, 1] = coefficient
    jacobian[5, 2] = coefficient
    jacobian[4, 4] = -1
    jacobian[5, 5] = -1
    jacobian[6, 3] = coefficient
    jacobian[6, 6] = -1

    jacobian[7, 1] = frequency
    jacobian[8, 2] = frequency
    jacobian[7:8, 3] .= -r_global .* y[6] .+ frequency .* d_global
    jacobian[7, 4] = 1
    jacobian[8, 5] = 1
    jacobian[7:8, 6] .= d_global
    return jacobian
end

function historical_first_order_stabilized_residual!(equations, t, y, dy,
                                                       p::BaumgarteParameters)
    implicit_first_order_stabilized_pendulum!(equations, dy, y, p, t)
    return nothing
end

function historical_first_order_stabilized_jacobian!(jacobian, t, y, dy,
                                                       coefficient,
                                                       p::BaumgarteParameters)
    jacobian .= implicit_first_order_stabilized_jacobian(
        t, y, dy, coefficient, p)
    return nothing
end

"""Construct initial data, optionally with prescribed constraint errors."""
function baumgarte_initial_conditions(theta, omega,
        p::BaumgarteParameters;
        position_error = [0.0, 0.0], velocity_error = [0.0, 0.0])
    mechanical = p.mechanical
    frequency = correction_frequency(p)
    r_global, d_global = marker_vectors(theta, mechanical)
    position = mechanical.pin - r_global + position_error
    velocity = -d_global .* omega + velocity_error

    K = zeros(5, 5)
    rhs = zeros(5)
    K[1, 1] = mechanical.mass
    K[2, 2] = mechanical.mass
    K[1, 4] = -1
    K[2, 5] = -1
    rhs[1:2] .= mechanical.mass .* mechanical.gravity
    K[3, 3] = mechanical.inertia
    K[3, 4:5] .= -d_global
    K[4, 1] = 1
    K[5, 2] = 1
    K[4:5, 3] .= d_global
    rhs[4:5] .= r_global .* omega^2 .-
        2frequency .* velocity_error .- frequency^2 .* position_error
    solution = K \ rhs

    y0 = [position; theta; velocity; omega; solution[4:5]]
    dy0 = [velocity; omega; solution[1:2]; solution[3]; 0.0; 0.0]
    return y0, dy0
end

"""Consistent initial data for first-order stabilized constraint errors."""
function first_order_stabilized_initial_conditions(theta, omega,
        p::BaumgarteParameters; position_error = [0.0, 0.0])
    mechanical = p.mechanical
    frequency = correction_frequency(p)
    r_global, d_global = marker_vectors(theta, mechanical)
    position = mechanical.pin - r_global + position_error
    phi_dot = -frequency .* position_error
    velocity = -d_global .* omega + phi_dot

    K = zeros(5, 5)
    rhs = zeros(5)
    K[1, 1] = mechanical.mass
    K[2, 2] = mechanical.mass
    K[1, 4] = -1
    K[2, 5] = -1
    rhs[1:2] .= mechanical.mass .* mechanical.gravity
    K[3, 3] = mechanical.inertia
    K[3, 4:5] .= -d_global
    K[4, 1] = 1
    K[5, 2] = 1
    K[4:5, 3] .= d_global
    rhs[4:5] .= r_global .* omega^2 .- frequency .* phi_dot
    solution = K \ rhs

    y0 = [position; theta; velocity; omega; solution[4:5]]
    dy0 = [velocity; omega; solution[1:2]; solution[3]; 0.0; 0.0]
    return y0, dy0
end

function baumgarte_problem(; theta0 = deg2rad(45.0), omega0 = 0.0,
        tspan = (0.0, 5.0), correction_time = 0.1,
        position_error = [0.0, 0.0], velocity_error = [0.0, 0.0])
    p = BaumgarteParameters(PendulumParameters(), correction_time)
    y0, dy0 = baumgarte_initial_conditions(theta0, omega0, p;
        position_error, velocity_error)
    function_definition = DAEFunction(implicit_baumgarte_pendulum!;
        jac = implicit_baumgarte_jacobian!)
    problem = DAEProblem(function_definition, dy0, y0, tspan, p;
        differential_vars = [trues(6); falses(2)])
    return problem
end

function run_baumgarte_ida(; theta0 = deg2rad(45.0), omega0 = 0.0,
        tspan = (0.0, 5.0), correction_time = 0.1,
        position_error = [1.0e-3, -5.0e-4],
        velocity_error = [0.0, 0.0], abstol = 1.0e-5,
        reaction_abstol = 1.0e-3, reltol = 1.0e-5)
    problem = baumgarte_problem(; theta0, omega0, tspan, correction_time,
        position_error, velocity_error)
    tolerances = [fill(abstol, 6); fill(reaction_abstol, 2)]
    return solve(problem, Sundials.IDA(); abstol = tolerances, reltol,
        initializealg = CheckInit())
end

function run_baumgarte_dassl(; theta0 = deg2rad(45.0), omega0 = 0.0,
        tspan = (0.0, 5.0), correction_time = 0.1,
        position_error = [1.0e-3, -5.0e-4],
        velocity_error = [0.0, 0.0], abstol = 1.0e-5,
        reaction_abstol = 1.0e-3, reltol = 1.0e-5,
        dt = 1.0e-5, dtmax = 0.05)
    p = BaumgarteParameters(PendulumParameters(), correction_time)
    y0, dy0 = baumgarte_initial_conditions(theta0, omega0, p;
        position_error, velocity_error)
    equations = zeros(8)
    residual = (t, y, dy) -> begin
        implicit_baumgarte_pendulum!(equations, dy, y, p, t)
        copy(equations)
    end
    jacobian = (t, y, dy, coefficient) ->
        implicit_baumgarte_jacobian(t, y, dy, coefficient, p)
    times, values, derivatives = DASSL.dasslSolve(residual, y0,
        collect(tspan); abstol = [fill(abstol, 6); fill(reaction_abstol, 2)],
        reltol = fill(reltol, 8), initstep = dt, maxstep = dtmax,
        dy0, jacobian)
    return (t = times, u = values, du = derivatives,
        retcode = ReturnCode.Success)
end


"""Run the stabilized pendulum with the Julia DDASSL reconstruction."""
function run_baumgarte_historical_ddassl(;
        theta0 = deg2rad(45.0), omega0 = 0.0,
        tspan = (0.0, 5.0), correction_time = 0.1,
        position_error = [1.0e-3, -5.0e-4],
        velocity_error = [0.0, 0.0], atol = 1.0e-5,
        rtol = 1.0e-5, dt = 1.0e-5, dtmax = 0.05)
    p = BaumgarteParameters(PendulumParameters(), correction_time)
    y0, dy0 = baumgarte_initial_conditions(theta0, omega0, p;
        position_error, velocity_error)
    options = HistoricalDDASSL.DASSLOptions{Float64}(
        atol = atol,
        rtol = rtol,
        initial_step = dt,
        maximum_step = dtmax,
    )

    variable_levels = [0, 0, 0, 1, 1, 1, 2, 2]
    equation_levels = [2, 2, 2, 1, 1, 1, 2, 2]
    differential_vars = [trues(6); falses(2)]

    return HistoricalDDASSL.dassl(
        historical_baumgarte_residual!, y0, dy0, tspan;
        parameter = p,
        jacobian! = historical_baumgarte_jacobian!,
        options,
        variable_levels,
        equation_levels,
        differential_vars,
        deficit = 0,
    )
end

"""Run first-order velocity-displacement stabilization."""
function run_first_order_stabilized_historical_ddassl(;
        theta0 = deg2rad(45.0), omega0 = 0.0,
        tspan = (0.0, 5.0), correction_time = 0.1,
        position_error = [1.0e-3, -5.0e-4], atol = 1.0e-5,
        rtol = 1.0e-5, dt = 1.0e-5, dtmax = 0.05)
    p = BaumgarteParameters(PendulumParameters(), correction_time)
    y0, dy0 = first_order_stabilized_initial_conditions(
        theta0, omega0, p; position_error)
    options = HistoricalDDASSL.DASSLOptions{Float64}(
        atol = atol,
        rtol = rtol,
        initial_step = dt,
        maximum_step = dtmax,
    )

    variable_levels = [0, 0, 0, 1, 1, 1, 2, 2]
    equation_levels = [2, 2, 2, 1, 1, 1, 1, 1]
    differential_vars = [trues(6); falses(2)]

    return HistoricalDDASSL.dassl(
        historical_first_order_stabilized_residual!, y0, dy0, tspan;
        parameter = p,
        jacobian! = historical_first_order_stabilized_jacobian!,
        options,
        variable_levels,
        equation_levels,
        differential_vars,
        deficit = 1,
    )
end

function baumgarte_diagnostics(solution, p::BaumgarteParameters)
    position_vectors = [position_constraint(y[1:3], p.mechanical)
                        for y in solution.u]
    velocity_vectors = [velocity_constraint(y[1:6], p.mechanical)
                        for y in solution.u]
    acceleration_vectors = map(solution.u, solution.du) do y, dy
        r_global, d_global = marker_vectors(y[3], p.mechanical)
        dy[4:5] + d_global .* dy[6] - r_global .* y[6]^2
    end
    position_errors = norm.(position_vectors, Inf)
    velocity_errors = norm.(velocity_vectors, Inf)
    acceleration_errors = norm.(acceleration_vectors, Inf)
    initial_position = first(position_vectors)
    initial_velocity = first(velocity_vectors)
    frequency = correction_frequency(p)
    initial_time = first(solution.t)
    decay_differences = map(solution.t, position_vectors) do t, phi
        elapsed = t - initial_time
        expected = (initial_position .+
            (initial_velocity .+ frequency .* initial_position) .* elapsed) .*
            exp(-frequency * elapsed)
        norm(phi - expected, Inf)
    end
    equation_errors = map(solution.u, solution.du, solution.t) do y, dy, t
        equations = zeros(eltype(y), length(y))
        implicit_baumgarte_pendulum!(equations, dy, y, p, t)
        norm(equations, Inf)
    end
    return (
        accepted_time_points = length(solution.t),
        initial_position_error = first(position_errors),
        final_position_error = last(position_errors),
        maximum_position_error = maximum(position_errors),
        final_velocity_error = last(velocity_errors),
        maximum_velocity_error = maximum(velocity_errors),
        final_acceleration_error = last(acceleration_errors),
        maximum_acceleration_error = maximum(acceleration_errors),
        maximum_critical_decay_difference = maximum(decay_differences),
        maximum_implicit_equation_error = maximum(equation_errors),
    )
end

function historical_baumgarte_solver_diagnostics(solution)
    nonzero_steps = abs.(solution.steps[2:end])
    return (
        accepted_steps = solution.stats.accepted_steps,
        rejected_steps = solution.stats.rejected_steps,
        residual_evaluations = solution.stats.residual_evaluations,
        jacobian_evaluations = solution.stats.jacobian_evaluations,
        factorizations = solution.stats.factorizations,
        newton_iterations = solution.stats.newton_iterations,
        maximum_order = maximum(solution.orders),
        minimum_step = minimum(nonzero_steps),
        maximum_step = maximum(nonzero_steps),
    )
end

function first_order_stabilization_diagnostics(solution,
                                                p::BaumgarteParameters)
    position_vectors = [position_constraint(y[1:3], p.mechanical)
                        for y in solution.u]
    velocity_vectors = [velocity_constraint(y[1:6], p.mechanical)
                        for y in solution.u]
    initial_position = first(position_vectors)
    frequency = correction_frequency(p)
    initial_time = first(solution.t)
    decay_differences = map(solution.t, position_vectors) do t, phi
        expected = initial_position .* exp(-frequency * (t - initial_time))
        norm(phi - expected, Inf)
    end
    equation_errors = map(solution.u, solution.du, solution.t) do y, dy, t
        equations = zeros(eltype(y), length(y))
        implicit_first_order_stabilized_pendulum!(equations, dy, y, p, t)
        norm(equations, Inf)
    end
    return (
        accepted_time_points = length(solution.t),
        initial_position_error = norm(initial_position, Inf),
        final_position_error = norm(last(position_vectors), Inf),
        maximum_velocity_error = maximum(norm.(velocity_vectors, Inf)),
        final_velocity_error = norm(last(velocity_vectors), Inf),
        maximum_exponential_decay_difference = maximum(decay_differences),
        maximum_implicit_equation_error = maximum(equation_errors),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    for correction_time in (0.2, 0.1, 0.05)
        println("correction_time = ", correction_time)
        p = BaumgarteParameters(PendulumParameters(), correction_time)
        for (name, runner) in (
                               ("Julia DDASSL reconstruction",
                                run_baumgarte_historical_ddassl),
                               ("DASSL", run_baumgarte_dassl),
                               ("IDA", run_baumgarte_ida))
            solution = runner(; correction_time)
            println("  ", name)
            for (key, value) in pairs(baumgarte_diagnostics(solution, p))
                println("    ", key, ": ", value)
            end
            if solution isa HistoricalDDASSL.DASSLResult
                for (key, value) in pairs(
                        historical_baumgarte_solver_diagnostics(solution))
                    println("    ", key, ": ", value)
                end
            end
        end
    end
    println("First-order velocity-displacement stabilization")
    for correction_time in (0.2, 0.1, 0.05)
        p = BaumgarteParameters(PendulumParameters(), correction_time)
        solution = run_first_order_stabilized_historical_ddassl(;
            correction_time)
        println("  correction_time = ", correction_time)
        for (key, value) in pairs(
                first_order_stabilization_diagnostics(solution, p))
            println("    ", key, ": ", value)
        end
        for (key, value) in pairs(
                historical_baumgarte_solver_diagnostics(solution))
            println("    ", key, ": ", value)
        end
    end
end
