@mtkmodel MarconatoModel begin
    @structural_parameters begin
        vf_input = true
        τ_m_input = true
    end
    @components begin
        terminal=Terminal()
        # inputs
        if vf_input
            vf_in = RealInput(guess=1) # field voltage input [pu]
        end
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
        X_d, [description="d-axis synchronous reactance"]
        X_q, [description="q-axis synchronous reactance"]
        X′_d, [description="d-axis transient reactance"]
        X′_q, [description="q-axis transient reactance"]
        X″_d, [description="d-axis subtransient reactance"]
        X″_q, [description="q-axis subtransient reactance"]
        X_ls, [description="stator leakage reactance"]
        T′_d0, [description="d-axis transient time constant"]
        T″_d0, [description="d-axis subtransient time constant"]
        T′_q0, [description="q-axis transient time constant"]
        T″_q0, [description="q-axis subtransient time constant"]
        T_AA, [description="d-axis additional leakage time constant"]
        H, [description="inertia constant"]
        D, [description="Damping constant"]
        # System and machine base
        S_b, [description="System power basis in MVA"]
        V_b, [description="System voltage basis in kV"]
        ω_b=1, [description="System base frequency in pu"]
        Ω_b, [description="Base synchronous frequency in rad/s"]
        Sn=S_b, [description="Machine power rating in MVA"]
        Vn=V_b, [description="Machine voltage rating in kV"]
        # input/parameter switches
        if !vf_input
            vf_set, [guess=1, description="field voltage"]
        end
        if !τ_m_input
            τ_m_set, [guess=1, description="mechanical torque"]
        end
    end
    @variables begin
        ψ_d(t), [description="d-axis flux linkage"]
        ψ_q(t), [description="q-axis flux linkage"]
        I_d(t), [description="d-axis current"]
        I_q(t), [description="q-axis current"]
        V_d(t), [description="d-axis voltage"]
        V_q(t), [description="q-axis voltage"]
        E′_d(t), [guess=0, description="transient voltage behind transient reactance in d-axis"]
        E′_q(t), [guess=1, description="transient voltage behind transient reactance in q-axis"]
        E″_d(t), [guess=0, description="subtransient voltage behind transient reactance in d-axis"]
        E″_q(t), [guess=1, description="subtransient voltage behind transient reactance in q-axis"]
        δ(t), [guess=0, description="rotor angle"]
        ω(t), [guess=1, description="rotor speed"]
        # observables
        v_mag(t), [description="terminal voltage [machine pu]"]
        v_arg(t), [description="Generator terminal angle"]
        P(t), [description="active power [machine pu]"]
        Q(t), [description="reactive power [machine pu]"]
        # inputs/parameters
        vf(t), [description="field voltage"]
        τ_m(t), [description="mechanical torque"]
        τ_e(t), [description="electrical torque"]
    end
    begin
        γ_d = T″_d0/T′_d0 * X″_d/X′_d * (X_d - X′_d)
        γ_q = T″_q0/T′_q0 * X″_q/X′_q * (X_q - X′_q)
        T_park(α) = [sin(α) cos(α); -cos(α) sin(α)]
    end
    @equations begin
        # Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_park(δ)*[V_d, V_q] * V_b/Vn
        # [terminal.i_r, terminal.i_i] .~ T_park(δ)*[I_d, I_q] * Ibase(S_b, V_b)/Ibase(Sn, Vn)
        # [V_d, V_q] .~ T_park(-δ)*[terminal.u_r, terminal.u_i] * Vn/V_b
        [I_d, I_q] .~ -T_park(-δ)*[terminal.i_r, terminal.i_i] * Ibase(Sn, Vn)/Ibase(S_b, V_b)

        0 ~ V_d + R_s * I_d + ψ_q
        0 ~ V_q + R_s * I_q - ψ_d
        
        0 ~ ψ_d + X″_d * I_d - E″_q #die beiden Gleichungen hat Sebastian nicht drin...
        0 ~ ψ_q + X″_q * I_q - E″_d

        Dt(E′_q) * T′_d0 ~ (-E′_q - (X_d - X′_d - γ_d)*I_d + (1 - T_AA/T′_d0) * vf) 
        Dt(E′_d) * T′_q0 ~ (-E′_d + (X_q - X′_q - γ_q)*I_q) 
        Dt(E″_q) * T″_d0 ~ (-E″_q + E′_q - (X′_d - X″_d - γ_d) * I_d + vf * T_AA/T′_d0)
        Dt(E″_d) * T″_q0 ~ (-E″_d + E′_d + (X′_q - X″_q - γ_q) * I_q)

        τ_e ~ ψ_d * I_q -  ψ_q * I_d #Sebastian: ((v_q + R_a * i_q) * i_q + (v_d + R_a * i_d) * i_d) / (ω + 1.0) 
        Dt(δ) ~ Ω_b * (ω - ω_b) # dθ ~ Ω * 2*pi * ω
        Dt(ω) * 2 * H ~ τ_m - τ_e - D * (ω - ω_b) #((P - D * ω)/(ω + 1.0) - τ_e) / (2*H)

        # inputs
        vf ~ vf_input ? vf_in.u : vf_set
        τ_m ~ τ_m_input ? τ_m_in.u : τ_m_set

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
