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

# compute vibrational partition function, given array of vibrational energies (with units of temperature)
function Z_vibr(T, E_vibr_array_K)
    return sum(exp.(-E_vibr_array_K ./ T))
end

# average some array over an equilibrium vibrational distribution, given array of vibrational energies (with units of temperature)
function avg_vibr_array(T, E_vibr_array_K, array_to_avg)
    Z_v = Z_vibr(T, E_vibr_array_K)

    return sum(array_to_avg .* exp.(-E_vibr_array_K ./ T)) / Z_v
end

# compute average vibrational energy, given array of vibrational energies (with units of temperature)
function e_vibr_from_array(m, T, E_vibr_array_K)
    return (k_B / m) * avg_vibr_array(T, E_vibr_array_K, E_vibr_array_K)
end

# compute non-scaled vibrational specific heat, cut-off harmonic oscillator, J/kg/K
function c_vibr_from_array(m, T, E_vibr_array_K)
    avg_e_sq = avg_vibr_array(T, E_vibr_array_K, E_vibr_array_K .^ 2)
    avg_e = avg_vibr_array(T, E_vibr_array_K, E_vibr_array_K)
    return (k_B / m) * (avg_e_sq - avg_e^2) / T^2
end

# generate array of vibrational energy (anharmonic spectrum) based on monotonicity criteria (e_i+1 > e_i)
function generate_e_vibr_arr_anharmonic_K_monotone(Θ, Θ_anh, E_diss)
    v_e_arr = []
    eold = -1.0
    i = 0
    while(true)
        enew = (i + 0.5) * Θ - (i + 0.5)^2 * Θ_anh
        if (enew > eold)
        eold = enew
        push!(v_e_arr, enew)
        i += 1
        else
        break
        end
    end
    return v_e_arr
end

function T_from_e_base(e_vibr_function, m, e_i, starting_T)
    target_f = T -> (5.0 / 2.0) * k_B * T / m .+ e_vibr_function(T) .- e_i
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
    points_shock = [1.288, 1.006, 2.15]
    cyl_radius = 1.0
    tanhparam1 = 1.4682760
    tanhparam2 = 1.9900576
    smult = 3
    sfactor = 1.02
    n_shock = 5
else
    points_shock = [1.2878, 1.0045, 2.148]
    cyl_radius = 1.0
    tanhparam1 = 1.0137246
    tanhparam2 = 1.5042155
    smult = 3
    sfactor = 1.02
    n_shock = 5
end

Θvibr = 2273.54  # O2
Θvibr_anh = 17.366  # O2
E_diss = 59364.8  # O2
mass = 5.3134e-26  # O2 mass
m_ref = mass  # O2 mass
R_specific = k_B / m_ref
rho_freestream = p_freestream / (R_specific * T_freestream)

L_ref = 0.045
p_ref = p_freestream
T_ref = T_freestream
rho_ref = rho_freestream
v_ref = sqrt(p_ref / rho_ref)

tmp_e_v_arr = generate_e_vibr_arr_anharmonic_K_monotone(Θvibr, Θvibr_anh, E_diss)

e_v_f = T -> e_vibr_from_array(mass, T, tmp_e_v_arr)
c_v_f = T -> c_vibr_from_array(mass, T, tmp_e_v_arr)
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