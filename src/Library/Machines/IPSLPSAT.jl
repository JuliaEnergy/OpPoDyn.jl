"""
  - Inputs V_f and output V_mag_out in system base not machien base
"""
@mtkmodel IPSLPSATOrder4 begin
    @components begin
        terminal = Terminal()
        # inputs
        P_m = RealInput() # mechanical power [pu]
        V_f = RealInput(guess=1) # field voltage input [pu]
        # outputs
        δ_out = RealOutput() # rotor angle
        ω_out = RealOutput() # rotor speed [pu]
        V_mag_out = RealOutput() # terminal voltage [pu]
        P_out = RealOutput() # active power [pu]
        Q_out = RealOutput() # reactive power [pu]
    end
    @parameters begin
        # base P
        R_a, [description="Armature resistance"]
        X′_d, [description="d-axis transient reactance"]
        M, [description="Mechanical starting time, 2H [Ws/VA]"]
        D, [description="Damping coefficient"]
        S_b, [description="System power basis in MVA"]
        V_b, [description="System voltage basis in kV"]
        ω_b, [description="System base frequency in rad/s"]
        S_n, [description="Machine power rating in MVA"]
        V_n, [description="Machine voltage rating in kV"]
        # 4th order params
        X_d=1.9, [description="d-axis synchronous reactance"]
        X_q=1.7, [description="q-axis synchronous reactance"]
        X′_q=0.5, [description="q-axis transient reactance"]
        T′_d0=8, [description="d-axis open circuit transient time constant"]
        T′_q0=0.8, [description="q-axis open circuit transient time constant"]
    end
    @variables begin
        # base vars
        V_arg(t), [description="Generator terminal angle"]
        V_d(t), [description="d-axis voltage"]
        V_q(t), [description="q-axis voltage"]
        I_d(t), [guess=0, description="d-axis current"]
        I_q(t), [guess=0, description="q-axis current"]
        P_e(t), [description="Electrical power transmitted through air gap"]

        δ(t), [guess=0, description="rotor angle"]
        ω(t), [guess=1, description="rotor speed [pu]"]
        V_mag(t), [description="terminal voltage [pu]"]
        P(t), [description="active power [pu]"]
        Q(t), [description="reactive power [pu]"]

        # 4th order vars
        V′_d(t), [guess=0, description="d-axis transient voltage"]
        V′_q(t), [guess=1, description="q-axis transient voltage"]
    end
    @equations begin
        # Park's transformations
        terminal.u_r ~ ( sin(δ) * V_d + cos(δ) *V_q) * V_b/V_n
        terminal.u_i ~ (-cos(δ) * V_d + sin(δ) * V_q) * V_b/V_n
        -terminal.i_r ~ (-sin(δ) * I_d - cos(δ) * I_q) * Ibase(S_b, V_b)/Ibase(S_n, V_n)
        -terminal.i_i ~ ( cos(δ) * I_d - sin(δ) * I_q) * Ibase(S_b, V_b)/Ibase(S_n, V_n)

        # observables
        V_mag ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        V_arg ~ atan(terminal.u_i, terminal.u_r)
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Q ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i

        #outputs
        P_out.u ~ P
        Q_out.u ~ Q
        V_mag_out.u ~ V_mag
        δ_out.u ~ δ
        ω_out.u ~ ω

        # swing equation
        P_e ~ (V_q + R_a*I_q)*I_q + (V_d + R_a*I_d)*I_d
        Dt(δ) ~ ω_b*(ω - 1)
        Dt(ω) ~ (P_m.u/S_b*S_n - P_e - D*(ω - 1))/M

        # internal transients
        Dt(V_′q) ~ ((-V′_q) - (X_d - X′_d)*I_d + V_f.u*V_b/V_n)/T′_d0;
        Dt(V′_d) ~ ifelse(abs(X_d - X′_q) < 1e-16,
                         ((-V′_d) + (X_q - X′_q)*I_q)/T′_q0,
                         (-V′_d)/T′_q0)
        V′_q ~ V_q + R_a*I_q + X′_d*I_d
        V′_d ~ V_d + R_a*I_d - X′_q*I_q
    end
end
