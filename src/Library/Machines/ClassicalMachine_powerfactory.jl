@mtkmodel ClassicalMachine_powerfactory begin
    #eher Modell aus Milano, das funktioniert so aber. Allerdings kann hier nicht τ einfach durch p_e ersetzt werden! Aber andere Vereinfachung (ω_b ~ 1 annzunehmen) ist auch notwendig!
    @structural_parameters begin
        P_m_input = true
    end
    @components begin
        terminal=Terminal()
        # inputs
        if P_m_input
            P_m_in = RealInput() # mechanical power [pu]
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
        H, [description="inertia constant in s related to Pgn"]
        D=0, [description="mechanical damping constant"]
        # System and machine base
        S_b, [description="System power basis in MVA"]
        V_b, [description="System voltage basis in kV"]
        ω_b, [description="System base frequency in rad/s"]
        #p=1, [description="Pole Pairs"]
        S_n=S_b, [description="Machine power rating in MVA"]
        V_n=V_b, [description="Machine voltage rating in kV"]
        # field voltage
        V_f_set, [guess=1, description="field voltage"] # = E´_q
        if !P_m_input
            P_m_set, [guess=1, description="mechanical power"]
        end
    end
    @variables begin
        δ(t), [guess=0, description="rotor angle"] #ableitung davon = omega
        ω(t), [guess=1, description="rotor speed"]
        I_d(t), [description="d-axis current"]
        I_q(t), [description="q-axis current"]
        V_d(t), [description="d-axis voltage"]
        V_q(t), [description="q-axis voltage"]
        P_m(t), [description="mechanical input power"]
        #P_e(t), [description="electrical output power"]
        #E´_q(t), [description="?"]
        #V_d_e(t), [description="internal d-axis voltage"]
        #V_q_e(t), [description="internal q-axis voltage"]
        # observables
        V_mag(t), [description="terminal voltage [machine pu]"]
        V_arg(t), [description="Generator terminal angle"]
        P(t), [description="active power [machine pu]"]
        Q(t), [description="reactive power [machine pu]"]
        #Ψ_d(t), [description="d-axis flux"]
        #Ψ_q(t), [description="q-axis flux"]
        # inputs/parameters
        #τ_m(t), [description="mechanical torque"]
        τ_e(t), [description="electrical torque"]
    end
    begin
        T_park(α) = [sin(α) cos(α); -cos(α) sin(α)]
    end
    @equations begin
        # Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_park(δ)*[V_d, V_q] * V_b/V_n
        #[terminal.i_r, terminal.i_i] .~ T_park(δ)*[I_d, I_q] * Ibase(S_b, V_b)/Ibase(Sn, Vn)
        #[V_d, V_q] .~ T_park(-δ)*[terminal.u_r, terminal.u_i] * Vn/V_b
        [I_d, I_q] .~ -T_park(-δ)*[terminal.i_r, terminal.i_i] * Ibase(S_n, V_n)/Ibase(S_b, V_b)

        # mechanical equations
        V_f_set - V_q ~ R_s * I_q +  X′_d * I_d
        V_d ~ X′_d * I_q - R_s * I_d
        2 * H * Dt(ω) ~ P_m / ω - τ_e - D * (ω - 1) #-ω_b -> funktioniert nur so!! (also wenn statt -ω_b -1 da steht!)
        Dt(δ) ~ ω_b * (ω - 1) #-ω_b > funktioniert nur so!! (also wenn statt -ω_b -1 da steht!) #/p #bis hier gleiche Gleichungen wie in Milano nur das Milano ω_b = 1 setzt (das macht aber keinen Unterschied)
        τ_e ~ (V_q + R_s*I_q)*I_q + (V_d + R_s*I_d)*I_d #nicht in openelectrical, sondern Milano


        # inputs
        P_m ~ P_m_input ? P_m_in.u : P_m_set # in openelectrical eigentlich P_m ~ (R_s + im * X′_d) ~ |V_t0|*|E_q0|*sin(δ_0 - θ_0)

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
