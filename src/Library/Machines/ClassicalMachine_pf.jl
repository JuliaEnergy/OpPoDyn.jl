@mtkmodel ClassicalMachine_pf begin
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
        δout = RealOutput() # rotor angle
        ωout = RealOutput() # rotor speed [pu]
        v_mag_out = RealOutput() # terminal voltage [pu]
        Pout = RealOutput() # active power [pu]
        Qout = RealOutput() # reactive power [pu]
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
        Sn=S_b, [description="Machine power rating in MVA"]
        Vn=V_b, [description="Machine voltage rating in kV"]
        # field voltage
        vf_set, [guess=1, description="field voltage"]
        if !τ_m_input
            τ_m_set, [guess=1, description="mechanical torque"]
        end
    end
    @variables begin
        δ(t), [guess=0, description="rotor angle"] #ableitung davon = omega
        ω(t), [guess=1, description="rotor speed"]
        I_d(t), [description="d-axis current"]
        I_q(t), [description="q-axis current"]
        V_d(t), [description="d-axis voltage"]
        V_q(t), [description="q-axis voltage"]
        #V_d_e(t), [description="internal d-axis voltage"] 
        #V_q_e(t), [description="internal q-axis voltage"] 
        # observables
        v_mag(t), [description="terminal voltage [machine pu]"]
        v_arg(t), [description="Generator terminal angle"]
        P(t), [description="active power [machine pu]"]
        Q(t), [description="reactive power [machine pu]"]
        #Ψ_d(t), [description="d-axis flux"] 
        #Ψ_q(t), [description="q-axis flux"] 
        # inputs/parameters
        #τ_m(t), [description="mechanical torque"]
        #τ_e(t), [description="electrical torque"]
        P_m(t), [description="mechanical input power"]
        #P_e(t), [description="electrical output power"]
        #E´(t), [description="voltage after transient reactance"]
    end
    begin
        T_park(α) = [sin(α) cos(α); -cos(α) sin(α)]
    end
    @equations begin
        # Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_park(δ)*[V_d, V_q] * V_b/Vn
        #[terminal.i_r, terminal.i_i] .~ T_park(δ)*[I_d, I_q] * Ibase(S_b, V_b)/Ibase(Sn, Vn)
        #[V_d, V_q] .~ T_park(-δ)*[terminal.u_r, terminal.u_i] * Vn/V_b
        [I_d, I_q] .~ -T_park(-δ)*[terminal.i_r, terminal.i_i] * Ibase(Sn, Vn)/Ibase(S_b, V_b)

        # mechanical swing equation
        2*H * Dt(ω) ~ P_m - P
        vf_set ~ Vn + I_d * im * X′_d
        P ~ Vn * vf_set * sin(δ)/X′_d
        Dt(δ) ~  ω

        # inputs
        #τ_m ~ τ_m_input ? τ_m_in.u : τ_m_set

        # observables
        v_mag ~ sqrt(V_d^2 + V_q^2)
        v_arg ~ atan(V_q, V_d)
        P ~ V_d*I_d + V_q*I_q
        Q ~ V_q*I_d - V_d*I_q
        
        #outputs
        Pout.u ~ P
        Qout.u ~ Q
        v_mag_out.u ~ v_mag
        δout.u ~ δ
        ωout.u ~ ω
    end
end