# A supersonic flow with gamma=1.4 around a cylinder
# using a StructuredMesh with an appropriate transformation
# so as to have shock-fitting and grid refinement near the shock

using OrdinaryDiffEq
using Trixi
using NLsolve
include("../src/compressible_euler_2d_vibrenergy.jl")


function mapping_full(xi_, eta_, n_ortho_points, cyl_radius, points_shock, n_shock, tanhparam1, tanhparam2,
    smult, sfactor)

    shock_pos = [(points_shock[1], 0.0), (points_shock[2], points_shock[2]), (0.0, points_shock[3])]  # 3 points that define shock

    # spline has form R[1] + c * eta_01^2 + d * eta_01^3, derivative w.r.t eta_01 is 0 at eta_01 = 0
    R = [sqrt(shock_pos[i][1]^2 + shock_pos[i][2]^2) for i in 1:3]  # 3 radii
    spline_matrix = [1.0 1.0; 0.25 0.125]  # find cubic spline coefficients
    spline_RHS = [R[3] - R[1], R[2] - R[1]]
    spline_cd = spline_matrix \ spline_RHS

    
    eta_01 = (eta_ + 1) / 2
    dxi = 2.0 / (n_ortho_points - 1)

    tanhparam = tanhparam1 + (tanhparam2 - tanhparam1) * (eta_01)  # linear interpolation between two tanhparams

    xi_shock = tanh(tanhparam * (2 - n_shock * dxi))
        
    xi_01 = tanh(tanhparam*(-xi_ + 1)) / xi_shock


    R_outer = R[1] + spline_cd[1] * eta_01^2 + spline_cd[2] * eta_01^3
    angle = -π/4 + eta_ * π/4

    multf = 1.0

    if (xi_ > -1.0 + smult * dxi)
        multf = 1.0
    else
        hbase = (smult * dxi)^2 
        multf = 1.0 + (xi_ - (-1.0 + smult * dxi))^2 * (sfactor - 1.0) / hbase
    end

    r = (cyl_radius + xi_01 * (R_outer - cyl_radius)) * multf

    # rounding to avoid nasty stuff like 1e-16 which makes line segment analysis in Paraview quite unpleasant
    return SVector(round(r * sin(angle); digits=8), round(r * cos(angle); digits=8))
end

@inline function initial_condition_supersonic_flow(x, t, equations::CompressibleEulerEquationsVibrEnergy2D)
    prim = SVector(rho_freestream / rho_ref, v1_freestream / v_ref, v2_freestream / v_ref, p_freestream / p_ref, T_freestream / T_ref)
    return prim2cons(prim, equations)
end

@inline function boundary_condition_supersonic_inflow(u_inner, normal_direction::AbstractVector, direction, x, t,
    surface_flux_function, equations::CompressibleEulerEquationsVibrEnergy2D)
    u_boundary = initial_condition_supersonic_flow(x, t, equations)

    flux = Trixi.flux(u_boundary, normal_direction, equations)

    return flux
end

@inline function boundary_condition_set(u_inner, normal_direction, direction, x, t,
    surface_flux_function, equations::CompressibleEulerEquationsVibrEnergy2D)
    # Calculate the boundary flux entirely from the internal solution state

    if direction == 1
        # inflow
        flux = boundary_condition_supersonic_inflow(u_inner, normal_direction, direction, x, t, surface_flux_function, equations)
    elseif direction == 2
        # cylinder surface
        flux = boundary_condition_slip_wall(u_inner, normal_direction, direction, x, t, surface_flux_function, equations)
    elseif direction == 3
        # axisymmetry
        flux = boundary_condition_slip_wall(u_inner, normal_direction, direction, x, t, surface_flux_function, equations)
    else 
        flux = Trixi.flux(u_inner, normal_direction, equations)  # supersonic outflow
    end

    return flux
end

function cons2prim_scaled(u, equations::CompressibleEulerEquationsVibrEnergy2D)
    # for output in non-scaled form, in SaveSolutionCallback: solution_variables=cons2prim_scaled
    rho, rho_v1, rho_v2, rho_e = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho

    e_internal = (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2)) / rho
    T = Trixi.temperature(e_internal, equations)
    p = rho * T

    return SVector(rho * rho_ref, v1 * v_ref, v2 * v_ref, p * p_ref, T * T_ref)
end

# compute non-scaled vibrational energy, infinite harmonic oscillator, J / kg
function e_vibr_iho(m, Θ, T)
    return (k_B / m) * Θ ./ ((exp.(Θ./T) .- 1.0))
end

# compute non-scaled vibrational specific heat, infinite harmonic oscillator, J/kg/K
function c_vibr_iho(m, Θ, T)
    return (k_B / m) * (Θ ./ T).^2 .* exp.(Θ ./ T) ./ (exp.(Θ ./ T) .- 1.0).^2
end

function T_from_e_base(e_vibr_function, m, e_i, starting_T)
    target_f = T -> (5.0 / 2.0) * k_B * T / m + e_vibr_function(T) .- e_i
    T_root = nlsolve(target_f, [starting_T])
    return T_root.zero[1]
end
  

Trixi.varnames(::typeof(cons2prim_scaled), ::CompressibleEulerEquationsVibrEnergy2D) = ("rho", "v1", "v2", "p", "T")

v1_freestream = 5956.0
v2_freestream = 0.0
p_freestream = 476.0
T_freestream = 901.0

tmax = 3.0  # 3.0-3.5 gives a converged solution

Nx = 30  # can be also 60
Ny = Nx
polydeg = 2


if Nx == 30
    points_shock = [1.2947, 1.0127, 2.163]
    cyl_radius = 1.0
    tanhparam1 = 1.4772730
    tanhparam2 = 1.9941675
    smult = 3
    sfactor = 1.02
    n_shock = 5
else
    points_shock = [1.294, 1.0115, 2.162]
    cyl_radius = 1.0
    tanhparam1 = 1.0217855
    tanhparam2 = 1.5083085
    smult = 3
    sfactor = 1.02
    n_shock = 5
end

Θvibr = 2273.54   # O2
mass = 5.3134e-26  # O2 mass
m_ref = mass  # O2 mass
R_specific = k_B / m_ref
rho_freestream = p_freestream / (R_specific * T_freestream)

L_ref = 0.045
p_ref = p_freestream
T_ref = T_freestream
rho_ref = rho_freestream
v_ref = sqrt(p_ref / rho_ref)

e_v_f = T -> e_vibr_iho(mass, Θvibr, T)
c_v_f = T -> c_vibr_iho(mass, Θvibr, T)
T_from_e = (e_i, T_i) -> T_from_e_base(e_v_f, mass, e_i, T_i)

e_ref = k_B * T_ref / mass
cv_ref = k_B / mass

cells_per_dimension = (Nx, Ny)

# needs to be in the Trixi namespace for VTK output
Trixi.mapping_full = mapping_full
Trixi.mapping_Nx = Nx
Trixi.cyl_radius = cyl_radius
Trixi.points_shock = points_shock
Trixi.tanhparam1 = tanhparam1
Trixi.tanhparam2 = tanhparam2
Trixi.smult = smult
Trixi.sfactor = sfactor
Trixi.n_shock = n_shock
mapping = (xi_, eta_) -> Trixi.mapping_full(xi_, eta_, Trixi.mapping_Nx, Trixi.cyl_radius,
                                            Trixi.points_shock, Trixi.n_shock, Trixi.tanhparam1, Trixi.tanhparam2,
                                            Trixi.smult, Trixi.sfactor)

mesh = StructuredMesh(cells_per_dimension, mapping)

equations = CompressibleEulerEquationsVibrEnergy2D(m_ref, e_v_f, c_v_f, T_from_e;
                                                   T_ref=T_ref, T_min=10.0, T_max=4.5e4, ΔT=1.0,
                                                   e_ref=e_ref)

initial_condition = initial_condition_supersonic_flow
boundary_conditions = boundary_condition_set

surface_flux = FluxLaxFriedrichs(max_abs_speed_naive_new) 
volume_flux  = FluxRotated(flux_oblapenko)

basis = LobattoLegendreBasis(polydeg)

indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max=0.5,
                                         alpha_min=0.001,
                                         alpha_smooth=true,
                                         variable=density_pressure)

volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux)

solver = DGSEM(polydeg=polydeg, surface_flux=surface_flux,
               volume_integral=volume_integral)


semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions=boundary_conditions)

tspan = (0.0, tmax)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

alive_callback = AliveCallback(alive_interval=100)
stepsize_callback = StepsizeCallback(cfl=0.7)

callbacks = CallbackSet(alive_callback,
                        stepsize_callback)


stage_limiter! = PositivityPreservingLimiterZhangShu(thresholds=(3.0e-4, 3.0e-4),
                                                     variables=(Trixi.density, pressure))
###############################################################################
# run the simulation


sol = solve(ode, SSPRK43(stage_limiter!);
            maxiters = 999999, ode_default_options()...,
            callback = callbacks);

summary_callback() # print the timer summary