using SciMLBase
using SparseArrays
using PracticalMechanicalSimulation.HistoricalDDASSL

@testset "DDASSL core behavior" begin
    residual! = (out, t, y, yprime, parameter) ->
        (out[1] = yprime[1] + parameter * y[1])
    jacobian! = (matrix, t, y, yprime, coefficient, parameter) ->
        (matrix[1, 1] = parameter + coefficient)
    accepted_step_times = Float64[]
    accepted_intervals = Tuple{Float64,Float64,Float64}[]
    accepted_step! = (time, y, yprime, order, step, stats) ->
        push!(accepted_step_times, time)
    accepted_interval! = function (left, right, nodes, values, order,
            step, stats)
        midpoint = (left + right) / 2
        value, _ = HistoricalDDASSL.interpolation_state(
            midpoint, nodes, values)
        push!(accepted_intervals, (left, right, value[1]))
    end
    solution = HistoricalDDASSL.dassl(residual!, [1.0], [-2.0], (0.0, 1.0);
        parameter = 2.0, jacobian!,
        options = HistoricalDDASSL.DASSLOptions{Float64}(
            atol = 1.0e-10, rtol = 1.0e-8,
            initial_step = 1.0e-5, maximum_step = 0.05),
        variable_levels = [1], equation_levels = [2],
        differential_vars = BitVector([true]),
        error_monitor = BitVector([true]), deficit = 0, accepted_step!,
        accepted_interval!)
    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution(1.0)[1] ≈ exp(-2.0) rtol = 2.0e-6
    @test solution.stats.accepted_steps > 0
    @test length(accepted_step_times) == solution.stats.accepted_steps
    @test accepted_step_times == solution.t[2:end]
    @test length(accepted_intervals) == solution.stats.accepted_steps
    @test last.(accepted_intervals) ≈ [
        solution((left + right) / 2)[1]
        for (left, right, _) in accepted_intervals]
    @test solution.stats.residual_evaluations > 0
    @test solution.stats.factorizations < solution.stats.accepted_steps
    @test solution.stats.jacobian_evaluations == solution.stats.factorizations
    @test solution.error_monitor.mask == Bool[true]
    @test length(solution.error_monitor.controlled_errors) == length(solution.t)
    @test length(solution.error_monitor.rms_errors) == length(solution.t)
    @test length(solution.error_monitor.maximum_errors) == length(solution.t)
    @test solution.error_monitor.maximum_error_indices[2:end] ==
        ones(Int, length(solution.t) - 1)

    diagnostic_times = Float64[]
    bounded = HistoricalDDASSL.dassl(residual!, [1.0], [-2.0], (0.0, 1.0);
        parameter = 2.0, jacobian!,
        options = HistoricalDDASSL.DASSLOptions{Float64}(
            atol = 1.0e-10, rtol = 1.0e-8,
            initial_step = 1.0e-5, maximum_step = 0.05),
        variable_levels = [1], equation_levels = [2],
        differential_vars = BitVector([true]),
        error_monitor = BitVector([true]), deficit = 0,
        save_everystep = false,
        accepted_diagnostics! = (time, y, yprime, order, step, stats,
            controlled, rms, maximum, index) ->
                push!(diagnostic_times, time))
    @test SciMLBase.successful_retcode(bounded.retcode)
    @test bounded.stats.accepted_steps > 2
    @test length(bounded.t) == 2
    @test length(bounded.y) == 2
    @test length(bounded.yprime) == 2
    @test length(bounded.error_monitor.maximum_errors) == 2
    @test last(bounded.t) == 1.0
    @test last(bounded.y)[1] ≈ exp(-2.0) rtol = 2.0e-6
    @test length(diagnostic_times) == bounded.stats.accepted_steps

    filtered_values = Float64[]
    filter_calls = Ref(0)
    filtered = HistoricalDDASSL.dassl(
        (out, t, y, yprime, parameter) -> (out[1] = yprime[1]),
        [0.0], [0.0], (0.0, 0.05);
        options = HistoricalDDASSL.DASSLOptions{Float64}(
            initial_step = 0.01, maximum_step = 0.02, maximum_order = 1),
        differential_vars = BitVector([true]),
        accepted_step_filter! =
            (time, y, yprime, order, step, stats) -> begin
                filter_calls[] += 1
                yprime[1] = 1
            end,
        accepted_step! = (time, y, yprime, order, step, stats) ->
            push!(filtered_values, yprime[1]))
    @test SciMLBase.successful_retcode(filtered.retcode)
    @test filter_calls[] == filtered.stats.accepted_steps
    @test filtered_values == first.(filtered.yprime[2:end])
    @test all(==(1), filtered.orders[2:end])

    # Dense output follows the corrected BDF history polynomial rather than a
    # separate cubic interpolant. A fourth-degree history distinguishes the
    # two methods while also checking its differentiated output.
    polynomial_times = collect(0.0:4.0)
    polynomial_values = [[time^4] for time in polynomial_times]
    polynomial_derivatives = [[4time^3] for time in polynomial_times]
    polynomial_solution = HistoricalDDASSL.DASSLResult(
        polynomial_times, polynomial_values, polynomial_derivatives,
        [1, 1, 2, 3, 4], [0.0; ones(4)],
        HistoricalDDASSL.DASSLStats(), SciMLBase.ReturnCode.Success,
        "test polynomial", nothing, nothing, nothing, BitVector([true]),
        nothing, 0,
        HistoricalDDASSL.DASSLEvent{Float64,Vector{Float64}}[])
    @test polynomial_solution(3.5)[1] ≈ 3.5^4 atol = 1.0e-12
    @test polynomial_solution(3.5, Val{1})[1] ≈ 4 * 3.5^3 atol = 1.0e-12
    @test polynomial_solution(3.0, Val{1})[1] == 4 * 3.0^3

    monitor_options = HistoricalDDASSL.DASSLOptions{Float64}(
        atol = 1.0, rtol = 0.0)
    monitor_summary = HistoricalDDASSL.monitored_error_summary(
        [0.0, 0.0], [-1.0, -1.0], monitor_options, [1.0, 0.1],
        Bool[true, true], 1)
    @test monitor_summary[1] ≈ sqrt((0.5^2 + 0.05^2) / 2)
    @test monitor_summary[2] ≈ 0.5
    @test monitor_summary[3] == 1

    failure_jacobian_calls = Ref(0)
    failure_jacobian! = function (matrix, t, y, yprime, coefficient, parameter)
        failure_jacobian_calls[] += 1
        matrix[1, 1] = failure_jacobian_calls[] <= 2 ?
            1.0e-6 : 1.0 + coefficient
    end
    sparse_prototype = sparse([1], [1], [0.0], 1, 1)
    recovered = HistoricalDDASSL.dassl(residual!, [1.0], [-1.0],
        (0.0, 0.2); parameter = 1.0, jacobian! = failure_jacobian!,
        jacobian_prototype = sparse_prototype,
        options = HistoricalDDASSL.DASSLOptions{Float64}(
            atol = 1.0e-8, rtol = 1.0e-6,
            initial_step = 0.1, maximum_step = 0.1,
            maximum_newton_iterations = 2,
            symbolic_reanalysis_failures = 2),
        variable_levels = [1], equation_levels = [2],
        differential_vars = BitVector([true]), deficit = 0)
    @test SciMLBase.successful_retcode(recovered.retcode)
    @test recovered(0.2)[1] ≈ exp(-0.2) rtol = 2.0e-4
    @test recovered.stats.corrector_failures == 2
    @test recovered.stats.symbolic_factorizations == 2
    @test recovered.stats.factorizations >
        recovered.stats.symbolic_factorizations

    changing_sign(t) = t < 0.15 ? 1.0 : -1.0
    changing_matrix! = (out, t, y, yprime, parameter) ->
        (out[1] = changing_sign(t) * (y[1] - cos(t)))
    changing_jacobian! = (matrix, t, y, yprime, coefficient, parameter) ->
        (matrix[1, 1] = changing_sign(t))
    refreshed = HistoricalDDASSL.dassl(changing_matrix!, [1.0], [0.0],
        (0.0, 0.3); jacobian! = changing_jacobian!,
        jacobian_prototype = sparse_prototype,
        options = HistoricalDDASSL.DASSLOptions{Float64}(
            atol = 1.0e-10, rtol = 1.0e-8,
            initial_step = 0.05, maximum_step = 0.05,
            maximum_order = 1,
            maximum_factorization_steps = 100_000),
        variable_levels = [0], equation_levels = [0],
        differential_vars = BitVector([false]),
        error_control = BitVector([true]), deficit = 0)
    @test SciMLBase.successful_retcode(refreshed.retcode)
    @test refreshed(0.3)[1] ≈ cos(0.3) atol = 1.0e-8
    @test refreshed.stats.factorizations == 2
    @test refreshed.stats.symbolic_factorizations == 1
    @test refreshed.stats.corrector_failures == 0

    singular_mode = Ref(false)
    singular_jacobian! =
        (matrix, t, y, yprime, coefficient, parameter) ->
            (matrix[1, 1] = singular_mode[] ? 1.0 + coefficient : 0.0)
    reconfiguration_reasons = Symbol[]
    reconfigure! = function (reason, t, y, yprime, error_control,
            differential_vars, parameter)
        push!(reconfiguration_reasons, reason)
        singular_mode[] = true
        (; jacobian_prototype = sparse([1], [1], [0.0], 1, 1),
           equation_levels = [1])
    end
    reconfigured = HistoricalDDASSL.dassl(residual!, [1.0], [-1.0],
        (0.0, 0.1); parameter = 1.0, jacobian! = singular_jacobian!,
        jacobian_prototype = sparse_prototype,
        variable_levels = [0], equation_levels = [1],
        differential_vars = BitVector([true]),
        error_control = BitVector([true]), reconfigure!)
    @test SciMLBase.successful_retcode(reconfigured.retcode)
    @test reconfigured(0.1)[1] ≈ exp(-0.1) rtol = 2.0e-6
    @test reconfiguration_reasons == [:singular_iteration_matrix]
    @test reconfigured.stats.state_reselections == 1

    algebraic! = (out, t, y, yprime, parameter) ->
        (out[1] = y[1] - cos(t))
    algebraic_jacobian! =
        (matrix, t, y, yprime, coefficient, parameter) ->
            (matrix[1, 1] = 1.0)
    algebraic = HistoricalDDASSL.dassl(algebraic!, [1.0], [0.0], (0.0, 1.0);
        jacobian! = algebraic_jacobian!,
        options = HistoricalDDASSL.DASSLOptions{Float64}(
            atol = 1.0e-10, rtol = 1.0e-8,
            initial_step = 1.0e-5, maximum_step = 0.05),
        variable_levels = [0], equation_levels = [0],
        differential_vars = BitVector([false]),
        error_control = BitVector([true]), deficit = 0)
    @test SciMLBase.successful_retcode(algebraic.retcode)
    @test algebraic(1.0)[1] ≈ cos(1.0) atol = 1.0e-8

    predictor_calls = Ref(0)
    coupled! = function (out, t, y, yprime, parameter)
        out[1] = yprime[1] + y[1]
        out[2] = y[2] - y[1]
    end
    coupled_jacobian! = function (matrix, t, y, yprime, coefficient,
            parameter)
        matrix[1, 1] = 1 + coefficient
        matrix[1, 2] = 0
        matrix[2, 1] = -1
        matrix[2, 2] = 1
    end
    predictor! = function (predicted, t, coefficient, history_derivative,
            parameter)
        predictor_calls[] += 1
        predicted[2] = predicted[1]
    end
    predicted = HistoricalDDASSL.dassl(coupled!, [1.0, 1.0], [-1.0, -1.0],
        (0.0, 1.0); jacobian! = coupled_jacobian!, predictor!,
        options = HistoricalDDASSL.DASSLOptions{Float64}(
            atol = 1.0e-10, rtol = 1.0e-8,
            initial_step = 1.0e-5, maximum_step = 0.05),
        variable_levels = [1, 0], equation_levels = [1, 0],
        differential_vars = BitVector([true, false]),
        error_control = BitVector([true, false]), deficit = 0)
    @test SciMLBase.successful_retcode(predicted.retcode)
    @test predictor_calls[] >= predicted.stats.accepted_steps
    @test predicted(1.0) ≈ fill(exp(-1), 2) rtol = 2.0e-6

    linear! = (out, t, y, yprime, parameter) ->
        (out[1] = yprime[1] - 1.0)
    linear_jacobian! =
        (matrix, t, y, yprime, coefficient, parameter) ->
            (matrix[1, 1] = coefficient)
    zero_surface! = (out, t, y, yprime, parameter) -> (out[1] = y[1])
    function linear_crossing(root_restart = :hard)
        HistoricalDDASSL.dassl(linear!, [-0.25], [1.0], (0.0, 1.0);
            jacobian! = linear_jacobian!, root! = zero_surface!,
            number_of_roots = 1, root_directions = [1], root_restart,
            options = HistoricalDDASSL.DASSLOptions{Float64}(
                atol = 1.0e-10, rtol = 1.0e-8,
                initial_step = 0.1, maximum_step = 0.2),
            variable_levels = [1], equation_levels = [2],
            differential_vars = BitVector([true]))
    end
    crossing = linear_crossing()
    @test SciMLBase.successful_retcode(crossing.retcode)
    @test length(crossing.events) == 1
    @test only(crossing.events).time ≈ 0.25 atol = 1.0e-12
    @test only(crossing.events).indices == [1]
    @test only(crossing.events).directions == [1]
    @test crossing.stats.events_found == 1
    @test crossing.stats.history_restarts == 1
    @test crossing.stats.root_evaluations > crossing.stats.accepted_steps
    soft_crossing = linear_crossing(:soft)
    no_restart_crossing = linear_crossing(:none)
    @test SciMLBase.successful_retcode(soft_crossing.retcode)
    @test SciMLBase.successful_retcode(no_restart_crossing.retcode)
    @test soft_crossing.stats.events_found == 1
    @test soft_crossing.stats.history_restarts == 1
    @test no_restart_crossing.stats.events_found == 1
    @test no_restart_crossing.stats.history_restarts == 0
    @test_throws ArgumentError linear_crossing(:invalid)

    two_surfaces! = function (out, t, y, yprime, parameter)
        out[1] = y[1]
        out[2] = y[1] - 0.2
    end
    two_crossings = HistoricalDDASSL.dassl(linear!, [-0.25], [1.0],
        (0.0, 0.6); jacobian! = linear_jacobian!, root! = two_surfaces!,
        number_of_roots = 2,
        options = HistoricalDDASSL.DASSLOptions{Float64}(
            atol = 1.0e-10, rtol = 1.0e-8,
            initial_step = 0.1, maximum_step = 0.3),
        variable_levels = [1], equation_levels = [2],
        differential_vars = BitVector([true]))
    @test getproperty.(two_crossings.events, :time) ≈ [0.25, 0.45] atol = 1e-12
    @test getproperty.(two_crossings.events, :indices) == [[1], [2]]

    filtered = HistoricalDDASSL.dassl(linear!, [-0.25], [1.0], (0.0, 0.5);
        jacobian! = linear_jacobian!, root! = zero_surface!,
        number_of_roots = 1, root_directions = [-1],
        options = HistoricalDDASSL.DASSLOptions{Float64}(
            initial_step = 0.1, maximum_step = 0.2),
        variable_levels = [1], equation_levels = [2],
        differential_vars = BitVector([true]))
    @test isempty(filtered.events)
end
