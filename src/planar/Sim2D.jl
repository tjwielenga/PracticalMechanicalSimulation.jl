"""
    Sim2D

Julia model-building interface for the planar mechanism simulator. Builders
create the same hierarchical model document used by TOML input, then pass it
through the ordinary planar loader and validation path.
"""
module Sim2D

import Base: run
import ..JuliaModelBuilder
using ..JuliaModelBuilder: document, element!, marker!, graphic!,
    set_properties!, analysis!, simulation!, state_selection!,
    initial_conditions!, parameters!, graphics!, variable, write_model,
    save_result
import ..PlanarModelIO: load_planar_model
import ..SimulationRunner: run_planar_model

export Model, ElementRef, VariableRef, document, element!, marker!, graphic!,
       set_properties!, analysis!, simulation!, state_selection!,
       initial_conditions!, parameters!, graphics!, variable, load_model,
       simulate, run, write_model, save_result,
       ground!, rigid_body!, floating_marker!, revolute!, inplane!, perp!,
       translational!, fixed!, span!, distance_coordinate!, coupler!,
       gear_pair!, rack_and_pinion!, pulley!, belt!, belt_span!, gravity!,
       applied_force!, applied_torque!, torsional_spring_damper!,
       spanning_force!, bushing!, plane_contact!, rotational_motion!,
       translational_motion!

const Model = JuliaModelBuilder.Model{:planar}
const ElementRef = JuliaModelBuilder.ElementRef{:planar}
const VariableRef = JuliaModelBuilder.VariableRef{:planar}

for (function_name, element_kind) in (
        (:ground!, "ground"),
        (:rigid_body!, "rigid_body"),
        (:floating_marker!, "floating_marker"),
        (:revolute!, "revolute"),
        (:inplane!, "inplane"),
        (:perp!, "perp"),
        (:translational!, "translational"),
        (:fixed!, "fixed"),
        (:span!, "span"),
        (:distance_coordinate!, "distance_coordinate"),
        (:coupler!, "coupler"),
        (:gear_pair!, "gear_pair"),
        (:rack_and_pinion!, "rack_and_pinion"),
        (:pulley!, "pulley"),
        (:belt!, "belt"),
        (:belt_span!, "belt_span"),
        (:gravity!, "gravity"),
        (:applied_force!, "applied_force"),
        (:applied_torque!, "applied_torque"),
        (:torsional_spring_damper!, "torsional_spring_damper"),
        (:spanning_force!, "spanning_force"),
        (:bushing!, "bushing"),
        (:plane_contact!, "plane_contact"),
        (:rotational_motion!, "rotational_motion"),
        (:translational_motion!, "translational_motion"))
    @eval begin
        function $(function_name)(model::Model, name; kwargs...)
            element!(model, $element_kind, name; kwargs...)
        end
    end
end

"""Validate and allocate a Julia-built model through the ordinary loader."""
function load_model(model::Model)
    load_planar_model(document(model);
        source_directory = model.source_directory)
end

function load_planar_model(model::Model;
        source_directory = model.source_directory)
    load_planar_model(document(model); source_directory)
end

"""Run a Julia-built model and return the ordinary Sim2D analysis result."""
function run_planar_model(model::Model; kwargs...)
    run_planar_model(load_model(model); kwargs...)
end

simulate(model::Model; kwargs...) = run_planar_model(model; kwargs...)
run(model::Model; kwargs...) = simulate(model; kwargs...)

end
