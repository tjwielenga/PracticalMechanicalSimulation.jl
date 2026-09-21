"""Dimension-independent scalar bristle friction laws."""
module FrictionLaws

export scalar_friction_coefficient, scalar_bristle_rate,
    scalar_bristle_force

scalar_friction_coefficient(friction, slip_squared) =
    friction.dynamic_coefficient +
    (friction.static_coefficient - friction.dynamic_coefficient) *
        exp(-slip_squared / friction.transition_speed^2)

function scalar_bristle_rate(friction, shear, slip, capacity)
    capacity <= 1.0e-12 && return -shear / friction.release_time
    slip - (friction.stiffness * abs(slip) / capacity) * shear
end

function scalar_bristle_force(friction, shear, slip, capacity)
    capacity <= 1.0e-12 && return zero(shear)
    rate = friction.stage == :static ? zero(shear) :
        scalar_bristle_rate(friction, shear, slip, capacity)
    trial = -friction.stiffness * shear - friction.damping * rate
    clamp(trial, -capacity, capacity)
end

end
