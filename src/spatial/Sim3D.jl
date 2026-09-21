"""
    Sim3D

Julia model-building interface for the spatial mechanism simulator. Builders
create the same hierarchical model document used by TOML and Lua input, then
pass it through the ordinary spatial loader and validation path.
"""
module Sim3D

import Base: run
import ..JuliaModelBuilder
using ..JuliaModelBuilder: document, element!, marker!, graphic!,
    set_properties!, analysis!, simulation!, state_selection!,
    initial_conditions!, parameters!, graphics!, variable, write_model,
    save_result
import ..SpatialModelIO: load_spatial_model
import ..SpatialSimulationRunner: run_spatial_model

export Model, ElementRef, VariableRef, document, element!, marker!, graphic!,
       set_properties!, analysis!, simulation!, state_selection!,
       initial_conditions!, parameters!, graphics!, variable, load_model,
       simulate, run, write_model, save_result,
       ground!, rigid_body!, gravity!, applied_force!, applied_torque!,
       spanning_force!, spherical!, perp!, inplane!, inline!, hinge!, orient!,
       revolute!, fixed!, rotational_motion!, translational_motion!,
       spanning_motion!, span!, directed_distance!, bushing!, plane_contact!,
       surface_friction!, revolute_friction!, translational_friction!,
       inplane_friction!, rolling_tire!, coupler!, gear_pair!,
       rack_and_pinion!, pulley!, belt!, belt_span!, equation_component!

const Model = JuliaModelBuilder.Model{:spatial}
const ElementRef = JuliaModelBuilder.ElementRef{:spatial}
const VariableRef = JuliaModelBuilder.VariableRef{:spatial}

for (function_name, element_kind) in (
        (:ground!, "ground"),
        (:rigid_body!, "rigid_body"),
        (:gravity!, "gravity"),
        (:applied_force!, "applied_force"),
        (:applied_torque!, "applied_torque"),
        (:spanning_force!, "spanning_force"),
        (:spherical!, "spherical"),
        (:perp!, "perp"),
        (:inplane!, "inplane"),
        (:inline!, "inline"),
        (:hinge!, "hinge"),
        (:orient!, "orient"),
        (:revolute!, "revolute"),
        (:fixed!, "fixed"),
        (:rotational_motion!, "rotational_motion"),
        (:translational_motion!, "translational_motion"),
        (:spanning_motion!, "spanning_motion"),
        (:span!, "span"),
        (:directed_distance!, "directed_distance"),
        (:bushing!, "bushing"),
        (:plane_contact!, "plane_contact"),
        (:surface_friction!, "surface_friction"),
        (:revolute_friction!, "revolute_friction"),
        (:translational_friction!, "translational_friction"),
        (:inplane_friction!, "inplane_friction"),
        (:rolling_tire!, "rolling_tire"),
        (:coupler!, "coupler"),
        (:gear_pair!, "gear_pair"),
        (:rack_and_pinion!, "rack_and_pinion"),
        (:pulley!, "pulley"),
        (:belt!, "belt"),
        (:belt_span!, "belt_span"),
        (:equation_component!, "equation_component"))
    @eval begin
        function $(function_name)(model::Model, name; kwargs...)
            element!(model, $element_kind, name; kwargs...)
        end
    end
end

"""Validate and allocate a Julia-built model through the ordinary loader."""
function load_model(model::Model)
    load_spatial_model(document(model);
        source_directory = model.source_directory,
        source_label = "Julia Sim3D model '$(model.document["model"]["name"])'")
end

function load_spatial_model(model::Model;
        source_directory = model.source_directory,
        source_label = "Julia Sim3D model '$(model.document["model"]["name"])'")
    load_spatial_model(document(model); source_directory, source_label)
end

"""Run a Julia-built model and return the ordinary Sim3D analysis result."""
function run_spatial_model(model::Model; kwargs...)
    run_spatial_model(load_model(model); kwargs...)
end

simulate(model::Model; kwargs...) = run_spatial_model(model; kwargs...)
run(model::Model; kwargs...) = simulate(model; kwargs...)

end
