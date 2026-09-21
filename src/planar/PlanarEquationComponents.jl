"""Planar adapter for the shared user-equation component implementation."""
module PlanarEquationComponents

using ..AutomaticAnalysis: EquationContribution
using ..EquationComponents
import ..PlanarComponentAssembly: executable_blocks, equation_contributions

export PlanarEquationComponent, planar_equation_registration,
    allocated_planar_equation_component, initialize_planar_equation_component!,
    planar_equation_state_rates!, set_planar_equation_stage!

const PlanarEquationComponent = UserEquationComponent
const planar_equation_registration = equation_component_registration
const allocated_planar_equation_component = allocated_equation_component
const initialize_planar_equation_component! = initialize_equation_component!
const planar_equation_state_rates! = equation_state_rates!
const set_planar_equation_stage! = set_equation_stage!

executable_blocks(component::PlanarEquationComponent) =
    equation_executable_blocks(component)
equation_contributions(::PlanarEquationComponent) = EquationContribution[]

end
