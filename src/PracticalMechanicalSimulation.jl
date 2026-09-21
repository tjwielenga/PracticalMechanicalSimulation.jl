"""
    PracticalMechanicalSimulation

Public entry point for the planar and spatial mechanism modelers. The package
reads hierarchical TOML or spatial Lua models, assembles their
component-local implicit equations into sparse canonical systems, runs the
requested analysis, and reads or writes portable `.simp` result files.
Spatial models may also be constructed directly in Julia with the exported
[`Sim3D`](@ref) builder module.

Most applications need only [`load_planar_model`](@ref),
[`run_planar_model`](@ref), [`load_spatial_model`](@ref),
[`run_spatial_model`](@ref), [`write_result`](@ref), and [`read_result`](@ref).
The included submodules remain accessible for element development and solver
experiments.
"""
module PracticalMechanicalSimulation

include("common/InputUnits.jl")
include("common/AssemblyExpansion.jl")
include("common/ScalarExpressions.jl")
include("common/AutomaticAnalysis.jl")
include("common/HistoricalDDASSL.jl")
include("common/ResultIO.jl")
include("common/SavedInitialConditions.jl")
include("common/ModalAnalysis.jl")
include("common/CSVExport.jl")
include("common/ModelExtraction.jl")

include("planar/PlanarAppliedForces.jl")
include("planar/PlanarDirectedDistances.jl")
include("planar/PlanarComponentAssembly.jl")
include("planar/PlanarModeling.jl")
include("planar/PlanarModelIO.jl")
include("planar/SimulationRunner.jl")
include("planar/CommandLine.jl")

include("spatial/SpatialComponentAssembly.jl")
include("spatial/SpatialModeling.jl")
include("spatial/SpatialDirectedDistances.jl")
include("spatial/SpatialConstraints.jl")
include("spatial/SpatialCoordinateCouplers.jl")
include("spatial/SpatialGearPairs.jl")
include("spatial/SpatialRackAndPinions.jl")
include("spatial/SpatialSpans.jl")
include("spatial/SpatialBelts.jl")
include("spatial/SpatialAppliedForces.jl")
include("spatial/SpatialBushings.jl")
include("spatial/SpatialPlaneContacts.jl")
include("spatial/SpatialFrictionForces.jl")
include("spatial/SpatialTires.jl")
include("spatial/SpatialEquationComponents.jl")
include("spatial/SpatialMotionGenerators.jl")
include("spatial/SpatialModelIO.jl")
include("spatial/SpatialSimulationRunner.jl")
include("spatial/SpatialCommandLine.jl")
include("spatial/Sim3D.jl")

using .PlanarModelIO: LoadedPlanarModel, PlanarModelMarker,
    load_planar_model, compile_time_expression
using .SimulationRunner: run_planar_model
using .ModalAnalysis: solve_modal_system
using .ResultIO: StoredBodyReference, StoredHealthSnapshot,
    StoredStateSelectionChange, StoredStaticSnapshot,
    StoredSimulationResult, IncrementalResultWriter,
    begin_incremental_result, record_incremental_event!,
    finish_incremental_result!, finalize_result!, write_result, read_result
using .CSVExport: export_result_csv, export_modal_csv, export_result_main
using .ModelExtraction: extract_model, extract_model_main
using .CommandLine: planar_model_main
using .SpatialModelIO: LoadedSpatialModel, load_spatial_model
using .SpatialSimulationRunner: run_spatial_model
using .SpatialCommandLine: spatial_model_main
import .Sim3D

export LoadedPlanarModel, PlanarModelMarker, load_planar_model,
       compile_time_expression,
       run_planar_model, solve_modal_system,
       StoredBodyReference, StoredHealthSnapshot, StoredStateSelectionChange,
       StoredStaticSnapshot, StoredSimulationResult, IncrementalResultWriter,
       begin_incremental_result, record_incremental_event!,
       finish_incremental_result!, finalize_result!, write_result, read_result,
       export_result_csv, export_modal_csv, export_result_main, extract_model,
       extract_model_main, planar_model_main,
       LoadedSpatialModel, load_spatial_model, run_spatial_model,
       spatial_model_main, Sim3D

end
