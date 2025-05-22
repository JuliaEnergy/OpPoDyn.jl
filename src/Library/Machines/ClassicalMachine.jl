@mtkmodel ClassicalMachine begin
    @structural_parameters begin
        τ_m_input = true
    end
    @components begin
        terminal=Terminal()
        # inputs
        if τ_m_input
            τ_m_in = RealInput() # mechanical torque [pu]
        end
        # outputs
        δ_out = RealOutput() # rotor angle
        ω_out = RealOutput() # rotor speed [pu]
        V_mag_out = RealOutput() # terminal voltage [pu]
        P_out = RealOutput() # active power [pu]
        Q_out = RealOutput() # reactive power [pu]
    end
    @parameters begin
        R_s, [description="stator resistance"]
        X′_d, [description="d-axis transient reactance"]
        H, [description="inertia constant"]
        D=0, [description="mechanical damping constant"]
        # System and machine base
        S_b, [description="System power basis in MVA"]
        V_b, [description="System voltage basis in kV"]
        ω_b, [description="System base frequency in rad/s"]
        S_n=S_b, [description="Machine power rating in MVA"]
        V_n=V_b, [description="Machine voltage rating in kV"]
        # field voltage
        V_f_set, [guess=1, description="field voltage"]
        if !τ_m_input
            τ_m_set, [guess=1, description="mechanical torque"]
        end
    end
    @variables begin
        δ(t), [guess=0, description="rotor angle"]
        ω(t), [guess=1, description="rotor speed"]
        I_d(t), [description="d-axis current"]
        I_q(t), [description="q-axis current"]
        V_d(t), [description="d-axis voltage"]
        V_q(t), [description="q-axis voltage"]
        # observables
        V_mag(t), [description="terminal voltage [machine pu]"]
        V_arg(t), [description="Generator terminal angle"]
        P(t), [description="active power [machine pu]"]
        Q(t), [description="reactive power [machine pu]"]
        # inputs/parameters
        τ_m(t), [description="mechanical torque"]
        τ_e(t), [description="electrical torque"]
    end
    begin
        T_park(α) = [sin(α) cos(α); -cos(α) sin(α)]
    end
    @equations begin
        # Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_park(δ)*[V_d, V_q] * V_b/V_n
        # [terminal.i_r, terminal.i_i] .~ T_park(δ)*[I_d, I_q] * Ibase(S_b, V_b)/Ibase(Sn, Vn)
        # [V_d, V_q] .~ T_park(-δ)*[terminal.u_r, terminal.u_i] * Vn/V_b
        [I_d, I_q] .~ -T_park(-δ)*[terminal.i_r, terminal.i_i] * Ibase(S_n, V_n)/Ibase(S_b, V_b)

        # mechanical swing equation
        Dt(δ) ~ ω_b*(ω - 1)
        2*H * Dt(ω) ~ τ_m / ω - τ_e - D*(ω-1)

        τ_e ~ (V_q + R_s*I_q)*I_q + (V_d + R_s*I_d)*I_d
        0 ~ V_q + R_s*I_q + X′_d*I_d - V_f_set
        0 ~ V_d + R_s*I_d - X′_d*I_q

        # inputs
        τ_m ~ τ_m_input ? τ_m_in.u : τ_m_set

        # observables
        V_mag ~ sqrt(V_d^2 + V_q^2)
        V_arg ~ atan(V_q, V_d)
        P ~ V_d*I_d + V_q*I_q
        Q ~ V_q*I_d - V_d*I_q

        #outputs
        P_out.u ~ P
        Q_out.u ~ Q
        V_mag_out.u ~ V_mag
        δ_out.u ~ δ
        ω_out.u ~ ω
    end
end
