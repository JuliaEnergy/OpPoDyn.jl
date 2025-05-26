@mtkmodel Swing begin
    @structural_parameters begin
        ω_ref_input=false;
        P_m_input=false;
    end
    @components begin
        terminal = Terminal()
        if ω_ref_input
            ω_ref = RealInput()
        end
        if P_m_input
            P_m = RealInput()
        end
    end
    @variables begin
        ω(t), [guess=0, description="Rotor frequency"]
        θ(t), [guess=0, description="Rotor angle"]
        P_e(t), [description="Electrical Power injected into the grid"]
    end
    @parameters begin
        M=0.005, [description="starting time, 2H [Ws/VA]"]
        D=0.0001, [description="Damping"]
        V_mag, [guess=1, description="Voltage magnitude"]
        if !ω_ref_input
            ω_ref=1, [description="Reference frequency"]
        end
        if !P_m_input
            P_m, [guess=1, description="Mechanical Power"]
        end
    end
    begin
        _ω_ref = ω_ref_input ? ω_ref.u : ω_ref
        _P_m = P_m_input ? P_m.u : P_m
    end
    @equations begin
        Dt(θ) ~ ω - _ω_ref
        M*Dt(ω) ~ _P_m - D*(ω - _ω_ref) - P_e

        P_e ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        terminal.u_r ~ V_mag*cos(θ)
        terminal.u_i ~ V_mag*sin(θ)
    end
end
