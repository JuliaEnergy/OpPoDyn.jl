@mtkmodel GovFixed begin
    @components begin
        τ_m = RealOutput()
    end
    @parameters begin
        τ_m_fixed, [guess=0, description="Fixed mechanical torque"]
    end
    @equations begin
        τ_m.u ~ τ_m_fixed
    end
end

_clamp(u, u_min, u_max) = max(min(u, u_max), u_min)

# from Milano P. 359
@mtkmodel TurbineGovTypeI begin
    @structural_parameters begin
        ω_ref_input=false
        p_ref_input=false
    end
    @components begin
        ω_meas = RealInput()
        if ω_ref_input
            ω_ref = RealInput()
        end
        if p_ref_input
            p_ref = RealInput()
        end
        τ_m = RealOutput()
    end
    @parameters begin
       if !ω_ref_input
            ω_ref=1, [description="Reference frequency"]
       end
       if !p_ref_input
            p_ref, [guess=1, description="Reference power [Machine PU]"]
       end
       p_min, [description="Minimum turbine output [Machine PU]"]
       p_max, [description="Maximum turbine output [Machine PU]"]
       R, [description="Govenor droop [Machine PU]"]
       # TODO: check Tc servo Ts governor
       T_c, [description="Servo time constant [s]"]
       T_s, [description="Govenor time constant [s]"]
       T_3, [description="Transient time constant 3 [s]"]
       T_4, [description="Transient time constant 4 [s]"]
       T_5, [description="Transient time constant 5 [s]"]
    end
    @variables begin
        p_droop(t), [description="P after droop (not limited)"]
        p_lim(t), [description="limited p"]
        x_g1(t)
        x_g2(t)
        x_g3(t)
    end
    begin
        _ω_ref = ω_ref_input ? ω_ref.u : ω_ref
        _p_ref = p_ref_input ? p_ref.u : p_ref
    end
    @equations begin
        p_droop ~ p_ref + 1/R * (_ω_ref - ω_meas)
        p_lim ~ _clamp(p_droop, p_min, p_max)
        T_s * Dt(x_g1) ~ p_lim - x_g1
        T_c * Dt(xg2) ~ (1-T_3/T_c)*x_g1 - x_g2
        T_5 * Dt(xg3) ~ (1-T_4/T_5)*(x_g2 + T_3/T_c*x_g1) - x_g3
        Dt(τ_m.u) ~ x_g3 + T_4/T_5*(x_g2 + T_3/T_c*x_g1)
    end
end

# from PSD.jl https://nrel-sienna.github.io/PowerSimulationsDynamics.jl/stable/component_models/turbine_gov/#TGOV1-[SteamTurbineGov1]
@mtkmodel TGOV1 begin
    @structural_parameters begin
        ω_ref_input=false
        p_ref_input=false
    end
    @components begin
        ω_meas = RealInput()
        if ω_ref_input
            ω_ref = RealInput()
        end
        if p_ref_input
            p_ref = RealInput()
        end
        τ_m = RealOutput()
    end
    @parameters begin
       if !ω_ref_input
            ω_ref=1, [description="Reference frequency"]
       end
       if !p_ref_input
            p_ref, [guess=1, description="Reference power [Machine PU]"]
       end
       V_min, [description="Valve min position"]
       V_max, [description="Valve max position"]
       R, [description="Govenor droop [Machine PU]"]
       T_1, [description="Transient time constant 1 [s]"]
       T_2, [description="Transient time constant 2 [s]"]
       T_3, [description="Transient time constant 3 [s]"]
       D, [description="Turbine Damping"]
    end
    @variables begin
        ref_sig(t), [description="Internal reference signal"]
        x_g1(t), [guess=1]
        # x_g1_sat(t)
        x_g2(t), [guess=0]
        Δω(t), [description="Speed deviation"]
    end
    begin
        _ω_ref = ω_ref_input ? ω_ref.u : ω_ref
        _p_ref = p_ref_input ? p_ref.u : p_ref
    end
    @equations begin
        # implementation after blockdiagram
        # https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Governor%20TGOV1%20and%20TGOV1D.htm
        # and milano

        # TODO: GOV: Δω absolute or relative?
        # Δω ~ ω_meas.u - _ω_ref
        Δω ~ ω_meas.u / _ω_ref - 1

        ref_sig ~ 1/R*(_p_ref  - Δω)

        T_1 * Dt(x_g1) ~ ifelse(
            ((x_g1 > V_max) & (ref_sig > x_g1)) | ((x_g1 < V_min) & (ref_sig < x_g1)),
            0,
            ref_sig - x_g1)

        T_3 * Dt(x_g2) ~ x_g1 + T_2*Dt(x_g1) - x_g2

        # TODO: GOV: output power or torque?
        τ_m.u*ω_meas.u ~ x_g2  - D*Δω
    end
end
