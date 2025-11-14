using Oceananigans
using Oceananigans.Units: day
using Oceananigans.Grids: znode, Center, AbstractTopology, Flat, Bounded
using Oceananigans.BoundaryConditions: ImpenetrableBoundaryCondition, fill_halo_regions!
using Oceananigans.Fields: ZeroField, ZFaceField
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry
using Adapt
import Adapt: adapt_structure, adapt
import Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_drift_velocity

const c = Center()

struct SUPRA{FT, FD, W} <: AbstractBiogeochemistry
    maximum_plankton_growth_rate    :: FT 
    maximum_bacteria_growth_rate    :: FT        
    bacteria_yield                  :: FT   
    zooplankton_yield               :: FT 
    zooplankton_grazing_coefficient :: FT  
    linear_mortality_rate           :: FT        
    quadratic_mortality_rate        :: FT  
    Z_quadratic_mortality_rate      :: FT  
    nutrient_half_saturation        :: FT     
    detritus_half_saturation        :: FT  
    incident_PAR                    :: FD    
    PAR_half_saturation             :: FT          
    PAR_attenuation_scale           :: FT        
    detritus_vertical_velocity      :: W       
end

"""
    SUPRA(; grid,
                                        maximum_plankton_growth_rate = 1/day,
                                        maximum_bacteria_growth_rate = 1/day
                                        maximum_grazing_rate         = 3/day
                                        bacteria_yield               = 0.2
                                        zooplankton_yield            = 0.3
                                        fraction_of_particulate_export = 0.33,
                                        linear_remineralization_rate = 0.03/day,
                                        linear_mortality_rate        = 0.01/day,
                                        quadratic_mortality_rate     = 0.1/day,
                                        quadratic_mortality_rate_Z   = 1/day,
                                        nutrient_half_saturation     = 0.1,
                                        detritus_half_saturation     = 0.1,
                                        grazing_half_saturation      = 3.0,
                                        PAR_half_saturation          = 10.0,
                                        PAR_attenuation_scale        = 25.0,
                                        detritus_vertical_velocity   = -10/day,
                                        depth_offset = -10)

Return a six-tracer biogeochemistry model for the interaction of nutrients (N), phytoplankton (P), 
zooplankton(Z), bacteria (B), dissolved detritus (D1), and particulate detritus (D2).

Keyword Arguments
=================
* `grid` (required): An Oceananigans' grid.

* `maximum_plankton_growth_rate`: (s⁻¹) Growth rate of plankton `P` unlimited by the
                                    availability of nutrients and light. Default: 1/day.

* `maximum_bacteria_growth_rate`: (s⁻¹) Growth rate of plankton `B` unlimited by the
                                  availability of nutrients and light. Default = 0.5/day.

* `maximum_grazing_rate`: (s⁻¹) Maximum grazing rate of phytoplankton by zooplankton.

* `bacteria_yield`: Determines fractional nutrient production by bacteria production 
                    relative to consumption of detritus such that ``∂_t N / ∂_t D = 1 - y``,
                    where `y = bacteria_yield`. Default: 0.2.

* `linear_remineralization_rate`: (s⁻¹) Remineralization rate constant of detritus 'D', 
                                  assuming linear remineralization of 'D', while 
                                  implicitly modeling bacteria 'B'. Default = 0.3/day.

* `linear_mortality_rate`: (s⁻¹) Linear term of the mortality rate of both plankton and bacteria.

* `quadratic_mortality_rate`: (s⁻¹) Quadratic term of the mortality rate of both plankton and bacteria.

* `nutrient_half_saturation`: (mmol m⁻³) Half-saturation of nutrients for plankton production.

* `detritus_half_saturation`: (mmol m⁻³) Half-saturation of nutrients for bacteria production.
                              Default = 10.0 mmol m⁻³.

* `phytoplankton_half_saturation`: (mmol m⁻³) Half-saturation of phytoplankton for zooplankton production.

* `zooplankton_assimilation`: Fractional assimilation efficiency for zooplankton.

* `PAR_half_saturation`: (W m⁻²) Half-saturation of photosynthetically available radiation (PAR)
                         for plankton production.

* `PAR_attenuation_scale`: (m) Depth scale over which photosynthetically available radiation (PAR)
                            attenuates exponentially.

* `detritus_sinking_speed`: (m s⁻¹) Sinking velocity of particulate detritus.

Tracer names
============
* `N`: nutrients

* `S`: superorganism

* `D`: detritus - particulate

Biogeochemical functions
========================
* transitions for `N`, `S`, `D`

* `biogeochemical_drift_velocity` for `D`, modeling the sinking of detritus at
  a constant `detritus_sinking_speed`.
"""
function SUPRA(; grid,
                maximum_plankton_growth_rate    = 1/day, # Add reference for each parameter
                maximum_bacteria_growth_rate    = 1/day,
                bacteria_yield                  = 0.2, 
                zooplankton_yield               = 0.3,
                zooplankton_grazing_coefficient = 0.5/day,
                linear_mortality_rate           = 0.01/day, # m³/mmol/day
                quadratic_mortality_rate        = 0.1/day,   # m³/mmol/day (zooplankton quadratic mortality)
                Z_quadratic_mortality_rate     = 1.0/day, 
                nutrient_half_saturation        = 0.1,      # mmol m⁻³
                detritus_half_saturation        = 0.1,      # mmol m⁻³
                incident_PAR                    = 700.0, # W m⁻²
                PAR_half_saturation             = 10.0,     # W m⁻²
                PAR_attenuation_scale           = 25.0,     # m
                detritus_vertical_velocity      = -10/day)  # m s⁻¹ )   # m 

    if detritus_vertical_velocity isa Number
        w₀ = detritus_vertical_velocity
        no_penetration = ImpenetrableBoundaryCondition()
        bcs = FieldBoundaryConditions(grid, (Center, Center, Face),
                                      top=no_penetration, bottom=no_penetration)
        detritus_vertical_velocity = ZFaceField(grid, boundary_conditions = bcs)
        set!(detritus_vertical_velocity, w₀)
        fill_halo_regions!(detritus_vertical_velocity)
    end

    if incident_PAR isa Number
        surface_PAR = incident_PAR            
        incident_PAR = CenterField(grid)            
        set!(incident_PAR, surface_PAR)            
        fill_halo_regions!(incident_PAR)
    elseif incident_PAR isa Field
        fill_halo_regions!(incident_PAR)
    end
    FD = typeof(incident_PAR)

    FT = eltype(grid)

    return SUPRA(convert(FT, maximum_plankton_growth_rate),   
                convert(FT, maximum_bacteria_growth_rate),            
                convert(FT, bacteria_yield),    
                convert(FT, zooplankton_yield),
                convert(FT, zooplankton_grazing_coefficient),      
                convert(FT, linear_mortality_rate),          
                convert(FT, quadratic_mortality_rate),  
                convert(FT, Z_quadratic_mortality_rate),      
                convert(FT, nutrient_half_saturation),       
                convert(FT, detritus_half_saturation),  
                incident_PAR,       
                convert(FT, PAR_half_saturation),            
                convert(FT, PAR_attenuation_scale),          
                detritus_vertical_velocity)
end

# const SUPRA = SUPRA

Adapt.adapt_structure(to, bgc::SUPRA) = 
SUPRA(adapt(to, bgc.maximum_plankton_growth_rate),   
    adapt(to, bgc.maximum_bacteria_growth_rate),        
    adapt(to, bgc.bacteria_yield),      
    adapt(to, bgc.zooplankton_yield),
    adapt(to, bgc.zooplankton_grazing_coefficient),  
    adapt(to, bgc.linear_mortality_rate),          
    adapt(to, bgc.quadratic_mortality_rate),
    adapt(to, bgc.Z_quadratic_mortality_rate),        
    adapt(to, bgc.nutrient_half_saturation),       
    adapt(to, bgc.detritus_half_saturation),   
    adapt(to, bgc.incident_PAR),   
    adapt(to, bgc.PAR_half_saturation),            
    adapt(to, bgc.PAR_attenuation_scale),          
    adapt(to, bgc.detritus_vertical_velocity))

@inline required_biogeochemical_tracers(::SUPRA) = (:N, :S, :Z, :D)

@inline function biogeochemical_drift_velocity(bgc::SUPRA, ::Val{:D})
    u = ZeroField()
    v = ZeroField()
    w = bgc.detritus_vertical_velocity
    return (; u, v, w)
end

# A depth-dependent temperature curve from Zakem (2018)
# Temp = 12 .*exp.(z./ 150) .+ 12 .*exp.(z ./ 500) .+ 2

# Temperature modification to metabolic rates, following the Arrhenius equation
# @inline temp_fun(Temp) = 0.8 .* exp.(-4000 .*(1 ./ (Temp .+ 273.15) .- 1 ./ 293.15))

# @inline bacteria_production(μᵇ, kᴰ, y, D, S) = y * μᵇ * D / (D + kᴰ) * S 
# @inline phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, P) = (μᵖ * min(N / (N + kᴺ) , I / (I + kᴵ)) * S) 


@inline function phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, S, μᵇ, kᴰ, y_b, D)
    BP_raw = y_b * μᵇ * D / (D + kᴰ)
    PP_raw = μᵖ * N / (N + kᴺ) * I / (I + kᴵ) 
    # ZP_raw = S * sqrt(max(0, (y_z * g * (BP_raw + PP_raw)))) - BP_raw - PP_raw
    ω_p = PP_raw / (PP_raw + BP_raw + eps()) 
    return ω_p * PP_raw * S
end

@inline function bacteria_production(μᵖ, kᴺ, kᴵ, I, N, S, μᵇ, kᴰ, y_b, D)
    BP_raw = y_b * μᵇ * D / (D + kᴰ) 
    PP_raw = μᵖ * N / (N + kᴺ) * I / (I + kᴵ) 
    # ZP_raw = S * sqrt(max(0, (y_z * g * (BP_raw + PP_raw)))) - BP_raw - PP_raw
    ω_b = BP_raw / (PP_raw + BP_raw + eps())   # eps() avoids division by zero
    return ω_b * BP_raw * S
end

@inline function omega_squared(μᵖ, kᴺ, kᴵ, I, N, μᵇ, kᴰ, y_b, D)
    BP_raw = y_b * μᵇ * D / (D + kᴰ) 
    PP_raw = μᵖ * N / (N + kᴺ) * I / (I + kᴵ) 
    ω_p = PP_raw / (PP_raw + BP_raw + eps()) 
    ω_b = BP_raw / (PP_raw + BP_raw + eps())   # eps() avoids division by zero
    return ω_p^2 + ω_b^2
end

@inline linear_mortality(mlin, S) = mlin * S 
@inline quadratic_mortality(mq, S) = mq * S^2
# @inline detritus_remineralization(r, D) = r * D

@inline function (bgc::SUPRA)(i, j, k, grid, ::Val{:N}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    μᵇ = bgc.maximum_bacteria_growth_rate
    kᴰ = bgc.detritus_half_saturation
    kᴺ = bgc.nutrient_half_saturation
    I₀ = bgc.incident_PAR
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    y_b = bgc.bacteria_yield
    y_z = bgc.zooplankton_yield
    g = bgc.zooplankton_grazing_coefficient

    # incoming shortwave
    z = znode(i, j, k, grid, c, c, c)
    I = I₀[i,j,k] * exp(z / λ)

    N = @inbounds fields.N[i, j, k]
    S = @inbounds fields.S[i, j, k]
    Z = @inbounds fields.Z[i, j, k]
    D = @inbounds fields.D[i, j, k] 

    # ω_z = zooplankton_omega(μᵖ, kᴺ, kᴵ, I, N, S, μᵇ, kᴰ, y_b, y_z, g, D)
    
    return (- phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, S, μᵇ, kᴰ, y_b, D) 
            + bacteria_production(μᵖ, kᴺ, kᴵ, I, N, S, μᵇ, kᴰ, y_b, D) * (1/y_b - 1) 
            + (1-y_z) * g * S * Z)
end

@inline function (bgc::SUPRA)(i, j, k, grid, ::Val{:S}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    kᴺ = bgc.nutrient_half_saturation
    I₀ = bgc.incident_PAR
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    mlin = bgc.linear_mortality_rate
    mq = bgc.quadratic_mortality_rate
    μᵇ = bgc.maximum_bacteria_growth_rate
    kᴰ = bgc.detritus_half_saturation
    y_b = bgc.bacteria_yield
    y_z = bgc.zooplankton_yield
    g = bgc.zooplankton_grazing_coefficient

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    # incoming shortwave
    I = I₀[i,j,k] * exp(z / λ)

    N = @inbounds fields.N[i, j, k]
    S = @inbounds fields.S[i, j, k]
    Z = @inbounds fields.Z[i, j, k]
    D = @inbounds fields.D[i, j, k] 

    return (phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, S, μᵇ, kᴰ, y_b, D) 
            + bacteria_production(μᵖ, kᴺ, kᴵ, I, N, S, μᵇ, kᴰ, y_b, D) 
            - g * S * Z
            - linear_mortality(mlin, S) 
            - quadratic_mortality(mq, S) * omega_squared(μᵖ, kᴺ, kᴵ, I, N, μᵇ, kᴰ, y_b, D))
end

@inline function (bgc::SUPRA)(i, j, k, grid, ::Val{:Z}, clock, fields)

    mlin = bgc.linear_mortality_rate
    mq_Z = bgc.Z_quadratic_mortality_rate
    y_z = bgc.zooplankton_yield
    g = bgc.zooplankton_grazing_coefficient

    S = @inbounds fields.S[i, j, k]
    Z = @inbounds fields.Z[i, j, k]

    return ( (y_z * g * S * Z)
            - linear_mortality(mlin, Z) - quadratic_mortality(mq_Z, Z))
end

@inline function (bgc::SUPRA)(i, j, k, grid, ::Val{:D}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    kᴺ = bgc.nutrient_half_saturation
    I₀ = bgc.incident_PAR
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    μᵇ = bgc.maximum_bacteria_growth_rate
    kᴰ = bgc.detritus_half_saturation
    y_b = bgc.bacteria_yield
    mlin = bgc.linear_mortality_rate
    mq = bgc.quadratic_mortality_rate
    mq_Z = bgc.Z_quadratic_mortality_rate

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    # incoming shortwave
    I = I₀[i,j,k] * exp(z / λ)

    N = @inbounds fields.N[i, j, k]
    D = @inbounds fields.D[i, j, k]
    S = @inbounds fields.S[i, j, k]
    Z = @inbounds fields.Z[i, j, k]

    return (linear_mortality(mlin, (S+Z)) 
            + quadratic_mortality(mq, S) * omega_squared(μᵖ, kᴺ, kᴵ, I, N, μᵇ, kᴰ, y_b, D)
            + quadratic_mortality(mq_Z, Z)
            - bacteria_production(μᵖ, kᴺ, kᴵ, I, N, S, μᵇ, kᴰ, y_b, D) / y_b)
end

