using PrecompileTools: @setup_workload, @compile_workload

# Lua assembly expansion and spatial-model construction contain a substantial
# amount of runtime dispatch. Compile one representative hierarchical model
# when the package is precompiled so that interactive loads in SimpView and at
# the Julia API do not pay that cost on their first use.
@setup_workload begin
    model_path = normpath(joinpath(@__DIR__, "..", "models", "spatial",
        "large-van.lua"))
    if isfile(model_path)
        @compile_workload begin
            load_spatial_model(model_path)
        end
    end
end
