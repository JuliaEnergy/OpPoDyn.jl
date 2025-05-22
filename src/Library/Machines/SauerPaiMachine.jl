@mtkmodel SauerPaiMachine begin
    @structural_parameters begin
        V_f_input = true
        τ_m_input = true
    end
    @components begin
        terminal=Terminal()
        # inputs
        if V_f_input
            V_f_in = RealInput(guess=1) # field voltage input [pu]
        end
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
        H, [description="inertia constant"]
        D=0, [description="direct shaft damping"]
        # System and machine base
        S_b, [description="System power basis in MVA"]
        V_b, [description="System voltage basis in kV"]
        ω_b, [description="System base frequency in rad/s"]
        S_n=S_b, [description="Machine power rating in MVA"]
        V_n=V_b, [description="Machine voltage rating in kV"]
        # input/parameter switches
        if !V_f_input
            V_f_set, [guess=1, bounds=(0,Inf), description="field voltage"]
        end
        if !τ_m_input
            τ_m_set, [guess=1, bounds=(0,Inf), description="mechanical torque"]
        end
    end
    @variables begin
        ψ_d(t), [description="d-axis flux linkage"]
        ψ_q(t), [description="q-axis flux linkage"]
        ψ″_d(t), [guess=1, description="flux linkage assosciated with X″_d"]
        ψ″_q(t), [guess=0, description="flux linkage assosciated with X″_q"]
        I_d(t), [guess=0, description="d-axis current"]
        I_q(t), [guess=0, description="q-axis current"]
        V_d(t), [guess=0, description="d-axis voltage"]
        V_q(t), [guess=1, description="q-axis voltage"]
        V′_d(t), [guess=1, description="transient voltage behind transient reactance in d-axis"]
        V′_q(t), [guess=0, description="transient voltage behind transient reactance in q-axis"]
        δ(t), [guess=0, description="rotor angle"]
        ω(t), [guess=1, description="rotor speed"]
        τ_e(t), [bounds=(0,Inf) ,description="electrical torque"]
        # observables
        V_mag(t), [description="terminal voltage [machine pu]"]
        V_arg(t), [description="Generator terminal angle"]
        P(t), [description="active power [machine pu]"]
        Q(t), [description="reactive power [machine pu]"]
        # inputs/parameters
        V_f(t), [bounds=(0,Inf), description="field voltage"]
        τ_m(t), [bounds=(0,Inf), description="mechanical torque"]
    end
    begin
        γ_d1 = (X″_d - X_ls)/(X′_d - X_ls)
        γ_q1 = (X″_q - X_ls)/(X′_q - X_ls)
        γ_d2 = (X′_d-X″_d)/(X′_d-X_ls)^2 # ~ (1 - γ_d1)/(X′_d - X_ls)
        γ_q2 = (X′_q-X″_q)/(X′_q-X_ls)^2 # ~ (1 - γ_q1)/(X′_q - X_ls)
        # This transformation seems identical to Milano and PSD models
        T_to_loc(α)  = [ sin(α) -cos(α);
                         cos(α)  sin(α)]
        T_to_glob(α) = [ sin(α)  cos(α);
                        -cos(α)  sin(α)]
    end
    @equations begin
        # Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_to_glob(δ)*[V_d, V_q] * V_n/V_b
        [I_d, I_q] .~ T_to_loc(δ)*[terminal.i_r, terminal.i_i] * Ibase(S_b, V_b)/Ibase(S_n, V_n)

        τ_e ~ ψ_d*I_q - ψ_q*I_d
        # for static ψ, this becomes which makes sense!
        # τ_e ~  (P + R_s*(I_d^2 + I_q^2))/ω

        Dt(δ) ~ ω_b*(ω - 1)
        2*H * Dt(ω) ~ τ_m  - τ_e - D*(ω - 1)

        # stator equations
        # 1/ω_b * Dt(ψ_d) ~ R_s*I_d + ω * ψ_q + V_d
        # 1/ω_b * Dt(ψ_q) ~ R_s*I_q - ω * ψ_d + V_q
        # static fomulation
        # 0 ~ R_s*I_d + ω * ψ_q + V_d
        # 0 ~ R_s*I_q - ω * ψ_d + V_q
        # static formualion in V_d, V_q which is the only free stuff
        V_d ~ -R_s*I_d - ω * ψ_q
        V_q ~ -R_s*I_q + ω * ψ_d

        T′_d0 * Dt(V′_q) ~ -V′_q - (X_d - X′_d)*(I_d - γ_d2*ψ″_d - (1-γ_d1)*I_d + γ_d2*V′_q) + V_f
        T′_q0 * Dt(V′_d) ~ -V′_d + (X_q - X′_q)*(I_q - γ_q2*ψ″_q - (1-γ_q1)*I_q - γ_q2*V′_d)
        T″_d0 * Dt(ψ″_d) ~ -ψ″_d + V′_q - (X′_d - X_ls)*I_d
        T″_q0 * Dt(ψ″_q) ~ -ψ″_q - V′_d - (X′_q - X_ls)*I_q

        # this constraint essentialy forces ψ_d and ψ_q
        ψ_d ~ -X″_d*I_d + γ_d1*V′_q + (1-γ_d1)*ψ″_d
        ψ_q ~ -X″_q*I_q - γ_q1*V′_d + (1-γ_q1)*ψ″_q
        # I_d ~ (-ψ_d + γ_d1*V′_q + (1-γ_d1)*ψ″_d)/X″_d
        # I_q ~ (-ψ_q - γ_q1*V′_d + (1-γ_q1)*ψ″_q)/X″_q
        # 0 ~ -X″_d*I_d + γ_d1*V′_q + (1-γ_d1)*ψ″_d - ψ_d
        # 0 ~ -X″_q*I_q - γ_q1*V′_d + (1-γ_q1)*ψ″_q - ψ_q

        # inputs
        V_f ~ V_f_input ? V_f_in.u : V_f_set
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


#=
MWE why whe need to have ψ force by constraint rather than free
in a nutshell: because the input i will force psi anyway, and if
and then the voltage depends on the derivative of I

@mtkmodel SauerPaiMachine2 begin
    @parameters begin
        I_d, [guess=0, description="d-axis current"]
        I_q, [guess=0, description="q-axis current"]
    end
    @variables begin
        ψ_d(t), [description="d-axis flux linkage"]
        ψ_q(t), [description="q-axis flux linkage"]
        V_d(t), [guess=0, description="d-axis voltage"]
        V_q(t), [guess=1, description="q-axis voltage"]
    end
    @equations begin
        Dt(ψ_d) ~ I_d + ψ_q + V_d
        Dt(ψ_q) ~ I_q - ψ_d + V_q
        0 ~ -I_d - ψ_d
        0 ~ -I_q - ψ_q
    end
end
VertexModel(Library.SauerPaiMachine2(name=:machine), [:I_d, :I_q], [:V_d, :V_q])
=#
