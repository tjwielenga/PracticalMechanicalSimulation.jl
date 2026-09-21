using Test

requested = isempty(ARGS) ? Set(["all"]) : Set(lowercase.(ARGS))
valid = Set(["all", "core", "planar", "spatial", "viewer", "ddassl"])
unknown = setdiff(requested, valid)
isempty(unknown) || error("unknown test suite: $(join(sort!(collect(unknown)), ", "))")
run_all = "all" in requested || "core" in requested

@testset "PracticalMechanicalSimulation" begin
    (run_all || "viewer" in requested) &&
        begin
            include("core/viewer_signal_browser_tests.jl")
            include("core/portable_viewer_document_tests.jl")
            include("core/simpview_server_tests.jl")
        end
    (run_all || "planar" in requested) &&
        begin
            include("core/planar_assembly_tests.jl")
            include("core/planar_julia_api_tests.jl")
            include("core/planar_equation_component_tests.jl")
            include("core/planar_friction_tests.jl")
            include("core/planar_stage_dependent_force_tests.jl")
            include("core/planar_curve_contact_tests.jl")
            include("core/model_program_tests.jl")
        end
    (run_all || "spatial" in requested) && begin
        include("core/spatial_julia_api_tests.jl")
        include("core/assembly_expansion_tests.jl")
        include("core/spatial_equation_component_tests.jl")
        include("core/spatial_surface_friction_tests.jl")
        include("core/spatial_revolute_friction_tests.jl")
        include("core/spatial_translational_friction_tests.jl")
        include("core/spatial_inplane_friction_tests.jl")
        include("core/spatial_tire_bristle_tests.jl")
        include("core/spatial_model_tests.jl")
    end
    (run_all || "ddassl" in requested) &&
        include("core/ddassl_tests.jl")
end
