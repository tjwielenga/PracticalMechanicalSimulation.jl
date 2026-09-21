"""Spatial adapter for the shared user-equation component implementation."""
module SpatialEquationComponents

using ..AutomaticAnalysis: EquationContribution
using ..EquationComponents
import ..SpatialComponentAssembly: executable_blocks, equation_contributions

export SpatialEquationComponent, spatial_equation_registration,
    allocated_spatial_equation_component, initialize_spatial_equation_component!,
    spatial_equation_state_rates!, set_spatial_equation_stage!

const SpatialEquationComponent = UserEquationComponent
const spatial_equation_registration = equation_component_registration
const allocated_spatial_equation_component = allocated_equation_component
const initialize_spatial_equation_component! = initialize_equation_component!
const spatial_equation_state_rates! = equation_state_rates!
const set_spatial_equation_stage! = set_equation_stage!

executable_blocks(component::SpatialEquationComponent) =
    equation_executable_blocks(component)
equation_contributions(::SpatialEquationComponent) = EquationContribution[]

end
