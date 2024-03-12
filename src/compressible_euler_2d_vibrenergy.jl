# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
using MuladdMacro
using Trixi
using StaticArrays
using LinearAlgebra

@muladd begin
    #! format: noindent
    
    @doc raw"""
    CompressibleEulerEquationsVibrEnergy2D(mass, e_vibr_function, c_vibr_function, T_from_e_function;
    T_ref=1.0, T_min=10.0, T_max=30.0e4, ΔT=1.0,
    e_ref=1.0, min_T_jump_rel=0.5)
    
    The compressible Euler equations
    ```math
    \frac{\partial}{\partial t}
    \begin{pmatrix}
    \rho \\ \rho v_1 \\ \rho v_2 \\ \rho e
    \end{pmatrix}
    +
    \frac{\partial}{\partial x}
    \begin{pmatrix}
     \rho v_1 \\ \rho v_1^2 + p \\ \rho v_1 v_2 \\ (\rho e +p) v_1
    \end{pmatrix}
    +
    \frac{\partial}{\partial y}
    \begin{pmatrix}
    \rho v_2 \\ \rho v_1 v_2 \\ \rho v_2^2 + p \\ (\rho e +p) v_2
    \end{pmatrix}
    =
    \begin{pmatrix}
    0 \\ 0 \\ 0 \\ 0
    \end{pmatrix}
    ```
    for an ideal gas with internal energy (per unit mass) given by
    ```math
    e_{int} = \frac{5}{2} \frac{kT}{m} + \mathrm{e\_vibr\_function}(T)
    ```
    in two space dimensions.
    Here, ``\rho`` is the density, ``v_1``, ``v_2`` the velocities, ``e`` the specific total energy **rather than** specific internal energy.
    """

    const k_B::Float64 = 1.380649e-23  # J / K

    struct CompressibleEulerEquationsVibrEnergy2D <:
           Trixi.AbstractCompressibleEulerEquations{2, 4}
        mass::Float64

        T_ref::Float64
        T_min::Float64
        T_max::Float64
        ΔT::Float64
        inv_ΔT::Float64
        min_T_jump::Float64

        e_arr::Vector{Float64}
        c_v_arr::Vector{Float64}
        R_specific::Float64  #  k_B / mass [J / kg / K] when not scaled

        e_ref::Float64
        c_v_ref::Float64

        # tabulated values used to find T(e) via linear interpolation
        e_min::Float64
        e_max::Float64
        Δe::Float64
        inv_Δe::Float64
        T_arr::Vector{Float64}

        # used to estimate \int c_v(tau) / tau d tau
        int_c_v_over_t_arr::Vector{Float64}
    
        # the e_vibr_function and c_vibr_function should compute usual dimensional quantities
        # scaling is governed by T_ref, e_ref, c_ref
        function CompressibleEulerEquationsVibrEnergy2D(mass, e_vibr_function, c_vibr_function, T_from_e_function;
                                                        T_ref=1.0, T_min=10.0, T_max=30.0e4, ΔT=1.0,
                                                        e_ref=1.0, min_T_jump_rel=0.5)
            n_range = trunc(Int, (T_max - T_min) / ΔT) + 1
            T_range = Vector(LinRange(T_min, T_max, n_range))
            
            c_v_ref = k_B / mass

            @assert abs((T_range[2] - T_range[1]) - ΔT) < 1e-3

            e_tot_from_T = T -> (5.0 / 2.0) * k_B * T / mass  .+ e_vibr_function(T)

            e_arr = zeros(n_range)
            c_v_arr = zeros(n_range)
            for i in 1:n_range
                e_arr[i] = e_tot_from_T(T_range[i])
            end

            c_v_from_T = T -> (5.0 / 2.0) * k_B / mass .+ c_vibr_function(T)


            for i in 1:n_range
                c_v_arr[i] = c_v_from_T(T_range[i])
            end

            e_min = e_arr[1]
            e_max = e_arr[n_range]

            # discretize energy range, will need to find corresponding temperature values at discretization points
            e_range = Vector(LinRange(e_min, e_max, n_range))
            Δe = e_range[2] - e_range[1]

            T_arr = zeros(n_range)

            # populate T(e) array
            for i in 1:n_range
                T_arr[i] = T_from_e_function(e_range[i], T_range[i])
            end

            int_c_v_over_t_arr = zeros(n_range)
            s_int_from_T = T -> c_v_from_T(T) / T
            # we integrate from T_min to T_max using Simpon's rule
            # int_{T_min}^{T_i} = int_{T_min}^{T_{i-1}} + int_{T_{i-1}}^{T_{i}}
            for i in 2:n_range
                T_a = T_range[i-1]
                T_b = T_range[i]

                int_c_v_over_t_arr[i] = int_c_v_over_t_arr[i-1]
                int_c_v_over_t_arr[i] += (ΔT / 6.0) * (s_int_from_T(T_a) + 4 * s_int_from_T(0.5 * (T_a + T_b)) + s_int_from_T(T_b))
            end

            T_min /= T_ref
            T_max /= T_ref
            ΔT /= T_ref
            T_arr ./= T_ref

            e_min /= e_ref
            e_max /= e_ref
            Δe /= e_ref
            e_arr ./= e_ref

            c_v_arr ./= c_v_ref
            int_c_v_over_t_arr ./= c_v_ref
            
            new(mass, T_ref, T_min, T_max, ΔT, 1.0 / ΔT, min_T_jump_rel * ΔT, e_arr, c_v_arr, (k_B / mass) / c_v_ref, e_ref, c_v_ref, e_min, e_max, Δe, 1.0 / Δe,
                T_arr, int_c_v_over_t_arr)
        end
    end
    
    function varnames(::typeof(cons2cons), ::CompressibleEulerEquationsVibrEnergy2D)
        ("rho", "rho_v1", "rho_v2", "rho_e")
    end

    function Trixi.varnames(::typeof(cons2cons), ::CompressibleEulerEquationsVibrEnergy2D)
        ("rho", "rho_v1", "rho_v2", "rho_e")
    end
    varnames(::typeof(cons2prim), ::CompressibleEulerEquationsVibrEnergy2D) = ("rho", "v1", "v2", "p", "T")
    

    # Set initial conditions at physical location `x` for time `t`
    """
        initial_condition_constant(x, t, equations::CompressibleEulerEquations2D)
    
    A constant initial condition to test free-stream preservation.
    """
    function initial_condition_constant(x, t, equations::CompressibleEulerEquationsVibrEnergy2D)
        rho = 1.0
        rho_v1 = 0.1
        rho_v2 = -0.2
        rho_e = 10.0
        return SVector(rho, rho_v1, rho_v2, rho_e)
    end
    
    """
        boundary_condition_slip_wall(u_inner, normal_direction, x, t, surface_flux_function,
                                     equations::CompressibleEulerEquationsVibrEnergy2D)
    
    Determine the boundary numerical surface flux for a slip wall condition.
    Imposes a zero normal velocity at the wall.
    Density is taken from the internal solution state and pressure is computed as an
    exact solution of a 1D Riemann problem. Further details about this boundary state
    are available in the paper:
    - J. J. W. van der Vegt and H. van der Ven (2002)
      Slip flow boundary conditions in discontinuous Galerkin discretizations of
      the Euler equations of gas dynamics
      [PDF](https://reports.nlr.nl/bitstream/handle/10921/692/TP-2002-300.pdf?sequence=1)
    
    Details about the 1D pressure Riemann solution can be found in Section 6.3.3 of the book
    - Eleuterio F. Toro (2009)
      Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction
      3rd edition
      [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)
    
    Should be used together with [`UnstructuredMesh2D`](@ref).
    """
    @inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
                                                  x, t,
                                                  surface_flux_function,
                                                  equations::CompressibleEulerEquationsVibrEnergy2D)
        norm_ = norm(normal_direction)
        # Normalize the vector without using `normalize` since we need to multiply by the `norm_` later
        normal = normal_direction / norm_
    
        # rotate the internal solution state
        u_local = Trixi.rotate_to_x(u_inner, normal, equations)
    
        # compute the primitive variables
        rho_local, v_normal, v_tangent, p_local, T_local = cons2prim(u_local, equations)
        gamma_local = get_gamma(T_local, equations)
    
        # Get the solution of the pressure Riemann problem
        # See Section 6.3.3 of
        # Eleuterio F. Toro (2009)
        # Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction
        # [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)
        # TODO: figure out???
        if v_normal <= 0.0
            sound_speed = sqrt(gamma_local * p_local / rho_local) # local sound speed
            p_star = p_local *
                     (1 + 0.5 * (gamma_local - 1) * v_normal / sound_speed)^(2 *
                                                                            gamma_local / (gamma_local - 1))
        else # v_normal > 0.0
            A = 2 / ((gamma_local + 1) * rho_local)
            B = p_local * (gamma_local - 1) / (gamma_local + 1)
            p_star = p_local +
                     0.5 * v_normal / A *
                     (v_normal + sqrt(v_normal^2 + 4 * A * (p_local + B)))
        end
    
        # For the slip wall we directly set the flux as the normal velocity is zero
        return SVector(zero(eltype(u_inner)),
                       p_star * normal[1],
                       p_star * normal[2],
                       zero(eltype(u_inner))) * norm_
    end
    
    """
        boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
                                     surface_flux_function, equations::CompressibleEulerEquations2D)
    
    Should be used together with [`TreeMesh`](@ref).
    """
    @inline function boundary_condition_slip_wall(u_inner, orientation,
                                                  direction, x, t,
                                                  surface_flux_function,
                                                  equations::CompressibleEulerEquationsVibrEnergy2D)
        # get the appropriate normal vector from the orientation
        if orientation == 1
            normal_direction = SVector(1, 0)
        else # orientation == 2
            normal_direction = SVector(0, 1)
        end
    
        # compute and return the flux using `boundary_condition_slip_wall` routine above
        return boundary_condition_slip_wall(u_inner, normal_direction, direction,
                                            x, t, surface_flux_function, equations)
    end
    
    """
        boundary_condition_slip_wall(u_inner, normal_direction, direction, x, t,
                                     surface_flux_function, equations::CompressibleEulerEquations2D)
    
    Should be used together with [`StructuredMesh`](@ref).
    """
    @inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
                                                  direction, x, t,
                                                  surface_flux_function,
                                                  equations::CompressibleEulerEquationsVibrEnergy2D)
        # flip sign of normal to make it outward pointing, then flip the sign of the normal flux back
        # to be inward pointing on the -x and -y sides due to the orientation convention used by StructuredMesh
        if isodd(direction)
            boundary_flux = -boundary_condition_slip_wall(u_inner, -normal_direction,
                                                          x, t, surface_flux_function,
                                                          equations)
        else
            boundary_flux = boundary_condition_slip_wall(u_inner, normal_direction,
                                                         x, t, surface_flux_function,
                                                         equations)
        end
    
        return boundary_flux
    end
    
    # Calculate 2D flux for a single point
    @inline function Trixi.flux(u, orientation::Integer, equations::CompressibleEulerEquationsVibrEnergy2D)
        rho, rho_v1, rho_v2, rho_e = u
        v1 = rho_v1 / rho
        v2 = rho_v2 / rho

        p = pressure(u, equations)
        if orientation == 1
            f1 = rho_v1
            f2 = rho_v1 * v1 + p
            f3 = rho_v1 * v2
            f4 = (rho_e + p) * v1
        else
            f1 = rho_v2
            f2 = rho_v2 * v1
            f3 = rho_v2 * v2 + p
            f4 = (rho_e + p) * v2
        end
        return SVector(f1, f2, f3, f4)
    end
    
    # Calculate 2D flux for a single point in the normal direction
    # Note, this directional vector is not normalized
    @inline function Trixi.flux(u, normal_direction::AbstractVector,
                          equations::CompressibleEulerEquationsVibrEnergy2D)
        rho_e = last(u)
        rho, v1, v2, p, T = cons2prim(u, equations)
    
        v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
        rho_v_normal = rho * v_normal
        f1 = rho_v_normal
        f2 = rho_v_normal * v1 + p * normal_direction[1]
        f3 = rho_v_normal * v2 + p * normal_direction[2]
        f4 = (rho_e + p) * v_normal
        return SVector(f1, f2, f3, f4)
    end

    @inline function flux_oblapenko(u_ll, u_rr, orientation::Integer,
                                    equations::CompressibleEulerEquationsVibrEnergy2D)
        # Unpack left and right state
        rho_ll, v1_ll, v2_ll, p_ll, T_ll = cons2prim(u_ll, equations)
        rho_rr, v1_rr, v2_rr, p_rr, T_rr = cons2prim(u_rr, equations)
        
        # Compute the necessary mean values
        rho_mean = Trixi.ln_mean(rho_ll, rho_rr)
        # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
        # in exact arithmetic since
        #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
        #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
        v1_avg = 0.5 * (v1_ll + v1_rr)
        v2_avg = 0.5 * (v2_ll + v2_rr)
        rho_avg = 0.5 * (rho_ll + rho_rr)
        inv_T_avg = 0.5 * (1.0 / T_ll + 1.0 / T_rr)
        T_geo = sqrt(T_ll * T_rr)

        e_int_ll = energy_internal_without_rho(u_ll, equations)
        e_int_rr = energy_internal_without_rho(u_rr, equations)

        e_internal_avg = 0.5 * (e_int_ll + e_int_rr)

        velocity_square_avg = 0.5 * (v1_ll^2 + v2_ll^2 + v1_rr^2 + v2_rr^2)

        T_jump = T_rr - T_ll

        if (abs(T_jump) < equations.min_T_jump)
            T_mid = 0.5 * (T_ll + T_rr)
            cv_Tast_over_Tast = c_v(T_mid, equations) / T_mid
            cv_T_astast = c_v(T_mid, equations)
        else
            int_cv_T_over_T_ll = entropy_c_v_integral(T_ll, equations)
            int_cv_T_over_T_rr = entropy_c_v_integral(T_rr, equations)
            cv_Tast_over_Tast = (int_cv_T_over_T_rr - int_cv_T_over_T_ll) / T_jump
            cv_T_astast = (e_int_rr - e_int_ll) / T_jump
        end

        # Calculate fluxes depending on orientation
        if orientation == 1
            f1 = rho_mean * v1_avg
            f2 = f1 * v1_avg + rho_avg / inv_T_avg
            f3 = f1 * v2_avg
            f4 = f1 * ((e_internal_avg - 0.5 * velocity_square_avg) + T_geo^2 * (cv_Tast_over_Tast - inv_T_avg * cv_T_astast)) +
                 f2 * v1_avg + f3 * v2_avg
        else
            f1 = rho_mean * v2_avg
            f2 = f1 * v1_avg
            f3 = f1 * v2_avg + rho_avg / inv_T_avg
            f4 = f1 * ((e_internal_avg - 0.5 * velocity_square_avg) + T_geo^2 * (cv_Tast_over_Tast - inv_T_avg * cv_T_astast)) +
                 f2 * v1_avg + f3 * v2_avg
        end
    
        return SVector(f1, f2, f3, f4)
    end
    
    # Calculate maximum wave speed for local Lax-Friedrichs-type dissipation as the
    # maximum velocity magnitude plus the maximum speed of sound
    @inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                         equations::CompressibleEulerEquationsVibrEnergy2D)
        rho_ll, v1_ll, v2_ll, p_ll, T_ll = cons2prim(u_ll, equations)
        rho_rr, v1_rr, v2_rr, p_rr, T_rr = cons2prim(u_rr, equations)
        gamma_ll = get_gamma(T_ll, equations)
        gamma_rr = get_gamma(T_rr, equations)
    
        # Get the velocity value in the appropriate direction
        if orientation == 1
            v_ll = v1_ll
            v_rr = v1_rr
        else # orientation == 2
            v_ll = v2_ll
            v_rr = v2_rr
        end
        # Calculate sound speeds
        c_ll = sqrt(gamma_ll * p_ll / rho_ll)
        c_rr = sqrt(gamma_rr * p_rr / rho_rr)
    
        return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr)
    end
    
    @inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                         equations::CompressibleEulerEquationsVibrEnergy2D)
        rho_ll, v1_ll, v2_ll, p_ll, T_ll = cons2prim(u_ll, equations)
        rho_rr, v1_rr, v2_rr, p_rr, T_rr = cons2prim(u_rr, equations)
        gamma_ll = get_gamma(T_ll, equations)
        gamma_rr = get_gamma(T_rr, equations)
    
        # Calculate normal velocities and sound speed
        # left
        v_ll = (v1_ll * normal_direction[1]
                +
                v2_ll * normal_direction[2])
        c_ll = sqrt(gamma_ll * p_ll / rho_ll)
        # right
        v_rr = (v1_rr * normal_direction[1]
                +
                v2_rr * normal_direction[2])
        c_rr = sqrt(gamma_rr * p_rr / rho_rr)
    
        return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
    end

    # Calculate maximum wave speed for local Lax-Friedrichs-type dissipation as the
    # maximum velocity magnitude plus the maximum speed of sound
    # better estimate than just naive
    @inline function max_abs_speed_naive_new(u_ll, u_rr, orientation::Integer,
                                             equations::CompressibleEulerEquationsVibrEnergy2D)
        rho_ll, v1_ll, v2_ll, p_ll, T_ll = cons2prim(u_ll, equations)
        rho_rr, v1_rr, v2_rr, p_rr, T_rr = cons2prim(u_rr, equations)
        gamma_ll = get_gamma(T_ll, equations)
        gamma_rr = get_gamma(T_rr, equations)
    
        # Get the velocity value in the appropriate direction
        if orientation == 1
            v_ll = v1_ll
            v_rr = v1_rr
        else # orientation == 2
            v_ll = v2_ll
            v_rr = v2_rr
        end
        # Calculate sound speeds
        c_ll = sqrt(gamma_ll * p_ll / rho_ll)
        c_rr = sqrt(gamma_rr * p_rr / rho_rr)
    
        λ_max = max(abs(v_ll) + c_ll, abs(v_rr) + c_rr)
        return λ_max
    end
    
    @inline function max_abs_speed_naive_new(u_ll, u_rr, normal_direction::AbstractVector,
                                             equations::CompressibleEulerEquationsVibrEnergy2D)
        rho_ll, v1_ll, v2_ll, p_ll, T_ll = cons2prim(u_ll, equations)
        rho_rr, v1_rr, v2_rr, p_rr, T_rr = cons2prim(u_rr, equations)
        gamma_ll = get_gamma(T_ll, equations)
        gamma_rr = get_gamma(T_rr, equations)
    
        # Calculate normal velocities and sound speed
        # left
        v_ll = (v1_ll * normal_direction[1]
                +
                v2_ll * normal_direction[2])
        c_ll = sqrt(gamma_ll * p_ll / rho_ll)
        # right
        v_rr = (v1_rr * normal_direction[1]
                +
                v2_rr * normal_direction[2])
        c_rr = sqrt(gamma_rr * p_rr / rho_rr)
    
        norm_norm = norm(normal_direction)
        return max(abs(v_ll) + c_ll * norm_norm, abs(v_rr) + c_rr * norm_norm)
    end
    
    # Called inside `FluxRotated` in `numerical_fluxes.jl` so the direction
    # has been normalized prior to this rotation of the state vector
    @inline function Trixi.rotate_to_x(u, normal_vector, equations::CompressibleEulerEquationsVibrEnergy2D)
        # cos and sin of the angle between the x-axis and the normalized normal_vector are
        # the normalized vector's x and y coordinates respectively (see unit circle).
        c = normal_vector[1]
        s = normal_vector[2]
    
        # Apply the 2D rotation matrix with normal and tangent directions of the form
        # [ 1    0    0   0;
        #   0   n_1  n_2  0;
        #   0   t_1  t_2  0;
        #   0    0    0   1 ]
        # where t_1 = -n_2 and t_2 = n_1
    
        return SVector(u[1],
                       c * u[2] + s * u[3],
                       -s * u[2] + c * u[3],
                       u[4])
    end
    
    # Called inside `FluxRotated` in `numerical_fluxes.jl` so the direction
    # has been normalized prior to this back-rotation of the state vector
    @inline function Trixi.rotate_from_x(u, normal_vector,
                                   equations::CompressibleEulerEquationsVibrEnergy2D)
        # cos and sin of the angle between the x-axis and the normalized normal_vector are
        # the normalized vector's x and y coordinates respectively (see unit circle).
        c = normal_vector[1]
        s = normal_vector[2]
    
        # Apply the 2D back-rotation matrix with normal and tangent directions of the form
        # [ 1    0    0   0;
        #   0   n_1  t_1  0;
        #   0   n_2  t_2  0;
        #   0    0    0   1 ]
        # where t_1 = -n_2 and t_2 = n_1
    
        return SVector(u[1],
                       c * u[2] - s * u[3],
                       s * u[2] + c * u[3],
                       u[4])
    end
    
    @inline function max_abs_speeds(u, equations::CompressibleEulerEquationsVibrEnergy2D)
        rho, v1, v2, p, gamma = cons2prim(u, equations)
        c = sqrt(gamma * p / rho)
    
        return abs(v1) + c, abs(v2) + c
    end
    
    # Convert conservative variables to primitive
    @inline function Trixi.cons2prim(u, equations::CompressibleEulerEquationsVibrEnergy2D)
        rho, rho_v1, rho_v2, rho_e = u

        v1 = rho_v1 / rho
        v2 = rho_v2 / rho

        e_internal = (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2)) / rho
        T = temperature(e_internal, equations)
        p = rho * T
    
        return SVector(rho, v1, v2, p, T)
    end

    # @inline function Trixi.get_gamma(T, equations::CompressibleEulerEquationsVibrEnergy2D)
    #     c_v_val = c_v(T, equations)
    #     return (c_v_val + 1.0) / c_v_val
    # end

    @inline function get_gamma(T, equations::CompressibleEulerEquationsVibrEnergy2D)
        c_v_val = c_v(T, equations)
        return (c_v_val + 1.0) / c_v_val
    end
    
    # Convert conservative variables to entropy
    @inline function Trixi.cons2entropy(u, equations::CompressibleEulerEquationsVibrEnergy2D)
        rho, v1, v2, p, T = cons2prim(u, equations)
        
        v_square = v1^2 + v2^2
        eint = energy_internal_without_rho(u, equations)

        s = entropy_thermodynamic(u, equations)
        # w1 = (equations.gamma - s) * equations.inv_gamma_minus_one - 0.5 * rho_p * v_square
        w1 = -s + (eint - 0.5 * v_square) / T
        w2 = v1 / T
        w3 = v2 / T
        w4 = -1.0 / T
    
        return SVector(w1, w2, w3, w4)
    end
    
    @inline function entropy2cons(w, equations::CompressibleEulerEquationsVibrEnergy2D)
        # See Hughes, Franca, Mallet (1986) A new finite element formulation for CFD
        # [DOI: 10.1016/0045-7825(86)90127-1](https://doi.org/10.1016/0045-7825(86)90127-1)
        # TODO
        return SVector(rho, rho_v1, rho_v2, rho_e)
    end
    
    # Convert primitive to conservative variables
    @inline function prim2cons(prim, equations::CompressibleEulerEquationsVibrEnergy2D)
        rho, v1, v2, p, T = prim
        rho_v1 = rho * v1
        rho_v2 = rho * v2
        rho_e = rho * energy(T, equations) + 0.5 * (rho_v1 * v1 + rho_v2 * v2)
        return SVector(rho, rho_v1, rho_v2, rho_e)
    end

    @inline function energy(T, equations::CompressibleEulerEquationsVibrEnergy2D)
        fracpos = (T - equations.T_min) * equations.inv_ΔT
        index_lower = floor(Int, fracpos)
        fracpos -= index_lower
        index_lower += 1
    
        return equations.e_arr[index_lower] * (1.0 - fracpos) + fracpos * equations.e_arr[index_lower + 1]
    end
    
    @inline function c_v(T, equations::CompressibleEulerEquationsVibrEnergy2D)
        fracpos = (T - equations.T_min) * equations.inv_ΔT
        index_lower = floor(Int, fracpos)
        fracpos -= index_lower
        index_lower += 1
    
        return equations.c_v_arr[index_lower] * (1.0 - fracpos) + fracpos * equations.c_v_arr[index_lower + 1]
    end
    
    @inline function entropy_c_v_integral(T, equations::CompressibleEulerEquationsVibrEnergy2D)
        fracpos = (T - equations.T_min) * equations.inv_ΔT
        index_lower = floor(Int, fracpos)
        fracpos -= index_lower
        index_lower += 1
    
        return equations.int_c_v_over_t_arr[index_lower] * (1.0 - fracpos) + fracpos * equations.int_c_v_over_t_arr[index_lower + 1]
    end
    
    @inline function Trixi.temperature(e, equations::CompressibleEulerEquationsVibrEnergy2D)
        if (e < equations.e_min)
            return equations.T_min
        end

        fracpos = (e - equations.e_min) * equations.inv_Δe
        index_lower = floor(Int, fracpos)
        fracpos -= index_lower
        index_lower += 1
        
        return equations.T_arr[index_lower] * (1.0 - fracpos) + fracpos * equations.T_arr[index_lower + 1]
    end

    @inline function temperature(e, equations::CompressibleEulerEquationsVibrEnergy2D)
        return Trixi.temperature(e, equations)
    end
    
    @inline function Trixi.density(u, equations::CompressibleEulerEquationsVibrEnergy2D)
        rho = u[1]
        return rho
    end
    
    @inline function Trixi.pressure(u, equations::CompressibleEulerEquationsVibrEnergy2D)
        rho, rho_v1, rho_v2, rho_e = u
        e_internal = (rho_e - 0.5 * (rho_v1^2 + rho_v2^2) / rho) / rho

        T = temperature(e_internal, equations)
        p = rho * T
        return p
    end
    
    @inline function Trixi.density_pressure(u, equations::CompressibleEulerEquationsVibrEnergy2D)
        rho, rho_v1, rho_v2, rho_e = u
    
        e_internal = (rho_e - 0.5 * (rho_v1^2 + rho_v2^2) / rho) / rho
        T = temperature(e_internal, equations)
        p = rho*T
        return rho * p
    end
    
    
    # Calculate thermodynamic entropy for a conservative state `cons`
    @inline function entropy_thermodynamic(cons, equations::CompressibleEulerEquationsVibrEnergy2D)
        # Pressure
        # p = (equations.gamma - 1) * (cons[4] - 1 / 2 * (cons[2]^2 + cons[3]^2) / cons[1])
    
        # Thermodynamic entropy
        # s = log(p) - equations.gamma * log(cons[1])
    
        e = energy_internal_without_rho(cons, equations)  # get e
        T = temperature(e, equations::CompressibleEulerEquationsVibrEnergy2D)
        s = entropy_c_v_integral(T, equations) - log(cons[1])
        return s
    end
    
    # Calculate mathematical entropy for a conservative state `cons`
    @inline function entropy_math(cons, equations::CompressibleEulerEquationsVibrEnergy2D)
        # Mathematical entropy
        S = -entropy_thermodynamic(cons, equations) * cons[1]
        return S
    end
    
    # Default entropy is the mathematical entropy
    @inline function entropy(cons, equations::CompressibleEulerEquationsVibrEnergy2D)
        entropy_math(cons, equations)
    end
    
    # Calculate total energy for a conservative state `cons`
    @inline energy_total(cons, ::CompressibleEulerEquationsVibrEnergy2D) = cons[4]
    
    # Calculate kinetic energy for a conservative state `cons`
    @inline function energy_kinetic(u, equations::CompressibleEulerEquationsVibrEnergy2D)
        rho, rho_v1, rho_v2, rho_e = u
        return (rho_v1^2 + rho_v2^2) / (2 * rho)
    end
    
    # Calculate internal energy for a conservative state `cons`
    @inline function energy_internal(cons, equations::CompressibleEulerEquationsVibrEnergy2D)
        # this returns rho e_internal [J/m^3]
        return energy_total(cons, equations) - energy_kinetic(cons, equations)
    end

    @inline function energy_internal_without_rho(cons, equations::CompressibleEulerEquationsVibrEnergy2D)
        # # this returns e_internal [J/kg]
        return return energy_internal(cons, equations) / cons[1]
    end
end # @muladd
    