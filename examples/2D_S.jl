# Script for initializing the 2D "AMOC" model with biogeochemistry 

########################################
########## Use Julia packages ##########
########################################
# using GLMakie
using CUDA
using Printf
using Statistics
using ClimaOceanBiogeochemistry: SUPRA
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: ZeroField, CenterField
using Oceananigans.BoundaryConditions: fill_halo_regions!

using Oceananigans.Models.HydrostaticFreeSurfaceModels:
                    HydrostaticFreeSurfaceModel,
                    PrescribedVelocityFields
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans: TendencyCallsite

########################################
######### Define model domain ##########
########################################

Ny = 500 
Nz = 200
Ly = 15000kilometers   # m
Lz = 4000           # m

arch = GPU()
# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(arch,
                       size = (Ny, Nz),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology=(Flat, Bounded, Bounded))

# Define streamfunction Ψ
deltaN = 500kilometers   # North downwelling width
deltaS = 3000kilometers  # South upwelling width
deltaZ = 2000     # Vertical asymmetry
Ψᵢ(y, z)  = - 7 * ((1 - exp(-y / deltaS)) * (1 - exp(-(Ly - y) / deltaN)) * 
            sinpi(z / Lz) * exp(z/deltaZ))

Ψ = Field{Center, Face, Face}(grid)
set!(Ψ, Ψᵢ)
fill_halo_regions!(Ψ, arch)

# Set velocity field from streamfunction
v = YFaceField(grid)
w = ZFaceField(grid)
v .= - ∂z(Ψ)
w .= + ∂y(Ψ)
fill_halo_regions!(v, arch)
fill_halo_regions!(w, arch)

# (I have to specify u to allow CheckPointer)
u = XFaceField(grid) 
fill_halo_regions!(u, arch)

########################################
############# Model setup ##############
########################################

# Vertical and horizontal diffusivity
kz(y,z,t) = 1e-4 + 5e-3 * (tanh((z+100)/20)+1) + 1e-2 * exp(-(z+4000)/50)
tracer_vertical_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), 
                                                κ=kz)
tracer_horizontal_closure = HorizontalScalarDiffusivity(κ=1e3)

# Set PAR as a function of latitude
incident_PAR =  Field{Nothing, Center, Center}(grid)
surface_PAR(y,z) = 100 + 600 * sinpi(y/Ly) 
set!(incident_PAR, surface_PAR)   
fill_halo_regions!(incident_PAR, arch)

# Model setup
model = HydrostaticFreeSurfaceModel(grid = grid,
                                    biogeochemistry = SUPRA(; grid,
                                                        incident_PAR = incident_PAR,),
                                    velocities = PrescribedVelocityFields(; u, v, w),
                                    tracers = (:N, :S, :Z, :D),
                                    tracer_advection = WENO(; order = 3),
                                    coriolis = nothing,
                                    buoyancy = nothing,
                                    closure = (tracer_vertical_closure, tracer_horizontal_closure))

# Initial conditions
set!(model, N=1.0, S=0.2, Z=0.1, D=0.01) # mol PO₄ m⁻³

# Model time
spinup_time = 365*2000days # 2000-year spin-up 
simulation = Simulation(model; Δt = 2hour, stop_time=spinup_time) 


# Print the progress 
# progress(sim) = @printf("Iteration: %d, time: %s, total(P): %.2e \n", 
#             iteration(sim), prettytime(sim),
#             sum(model.tracers.PO₄) + sum(model.tracers.POP) + sum(model.tracers.DOP))
# add_callback!(simulation, progress, IterationInterval(100))

# Save model results in the output file
# outputs = ( # v = model.velocities.v,
#             PO₄= model.tracers.PO₄,
#             Dremin = model.biogeochemistry.Dremin)

simulation.output_writers[:simple_output] =
        JLD2OutputWriter(model, model.tracers; 
                        schedule = TimeInterval(365*100days), 
                        filename = "2D_S",
                        overwrite_existing = false)

simulation.output_writers[:checkpointer] = Checkpointer(model,
            schedule = TimeInterval(compute_time),
            prefix = "2D_S_checkpoint",
            overwrite_existing = false)
        
run!(simulation, pickup = false)

#################################### Visualize ####################################
#=
# filepath = simulation.output_writers[:simple_output].filepath
filepath = "./2D_S.jld2"

N_timeseries = FieldTimeSeries(filepath, "N")
times = N_timeseries.times
xt, yt, zt = nodes(N_timeseries)
S_timeseries = FieldTimeSeries(filepath, "S")
Z_timeseries = FieldTimeSeries(filepath, "Z")
D_timeseries = FieldTimeSeries(filepath, "D")

n = Observable(1)
# title = @lift @sprintf("t = Year %d x 50", times[$n] / (50*365.25days)) 
title = @lift @sprintf("t = Year %d", times[$n] / 365days) 

# convert unit from mol/m³ to μM: 1e3*interior(...)
Nₙ = @lift interior(N_timeseries[$n], 1, :, :)
Sₙ = @lift interior(S_timeseries[$n], 1, :, :)
Zₙ = @lift interior(Z_timeseries[$n], 1, :, :)
Dₙ = @lift interior(D_timeseries[$n], 1, :, :)
# avg_Nₙ = @lift mean(interior(N_timeseries[$n], 1, :, :), dims=1) 

#################################################################################
fig = Figure(size=(1800, 500))

ax_N = Axis(fig[2, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "N", aspect = 1)
hm_N = heatmap!(ax_N, yt/1e3, zt, Nₙ; colorrange = (0,3), interpolate = false) 
Colorbar(fig[2, 2], hm_N; flipaxis = false)

ax_S = Axis(fig[2, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "S", aspect = 1)
hm_S = heatmap!(ax_S, yt/1e3, zt, Sₙ; colorrange = (0,0.5), interpolate = false) 
Colorbar(fig[2, 4], hm_S; flipaxis = false)

ax_Z = Axis(fig[2, 5]; xlabel = "y (km)", ylabel = "z (m)", title = "Z", aspect = 1)
hm_Z = heatmap!(ax_Z, yt/1e3, zt, Zₙ; colorrange = (0,0.5), interpolate = false) 
Colorbar(fig[2, 6], hm_Z; flipaxis = false)

ax_D = Axis(fig[2, 7]; xlabel = "y (km)", ylabel = "z (m)", title = "D", aspect = 1)
hm_D = heatmap!(ax_D, yt/1e3, zt, Dₙ; colorrange = (0,0.1),interpolate = false) 
Colorbar(fig[2, 8], hm_D; flipaxis = false)

# ax_avg_PO4 = Axis(fig[3, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", title = "Average [PO₄] (μM)", yaxisposition = :right)
# xlims!(ax_avg_PO4, 0, 3)
# PO4_prof = lines!(ax_avg_PO4, avg_PO4ₙ[][1, :], zt)

fig[1, 1:8] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "2D_S.mp4", frames, framerate=10) do i
    n[] = i
#     PO4_prof[1] = avg_PO4ₙ[][1, :]
end
nothing #hide
=#
