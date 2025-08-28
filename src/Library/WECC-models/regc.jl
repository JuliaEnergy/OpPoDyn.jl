#without High Voltage Reactive Current Managment and without Low Voltage Active Current Managment

@mtkmodel regc_a begin
    @components begin
        # inputs
        Vt_in = RealInput(guess=1)
        Iqcmd_in = RealInput(guess=1)
        Ipcmd_in = RealInput(guess=1)
        Q_gen0 = RealInput(guess=1) #woher?
        # outputs
        Iqout = RealOutput()
        Ipout = RealOutput()
    end
    @parameters begin
        I_qrmax, [description="Maximum rate-of-change of reactive current (pu/s)"]
        I_qrmin, [description="Minimum rate-of-change of reactive current (pu/s)"]
        T_g, [description="Inverter current regulator lag time constant (s)"]
        T_fltr, [description="Terminal voltage filter (for LVPL) time constant (s)"]
        Brkpt, [description="LVPL breakpoint (pu voltage)"]
        Zerox, [description="LVPL zero crossing (pu voltage)"]
        L_vpl1, [description="LVPL gain breakpoint (pu current on mbase / pu voltage)"]
        L_vplsw, [description="Enable (1) or disable (0) low voltage power logic"]
        rrpwr, [description="Active current up-ramp rate limit on voltage recovery (pu/s)"]
    end
    @variables begin
        I_qz(t), [description="I_q after inverter current regulator"]
        I_q(t), [description="I_q after inverter current regulator with rate limits"]
        I_pz(t), [description="I_p after Inverter current regulator"]
        I_p(t), [description="I_p after inverter current regulator with limits"]
        V(t), [description="V_t after filter"]
        I_lvpl(t), [description="current resulting from low voltage power logic"]
        L_vplsw_adjust(t), [description="switch low voltage power logic off if V > breakpoint"]
        I_p_target(t), [description="target I_p after LVPL limiting"]
    end
    @equations begin
        #q-phase current
        T_g * Dt(I_qz) ~ -Iqcmd_in.u - I_qz
        Dt(I_q) ~ ifelse(Q_gen0.u>0, uplimit(Dt(I_qz), I_qrmax), lowlimit(Dt(I_qz), I_qrmin))
        #p-phase current
        T_fltr * Dt(V) ~ Vt_in.u - V
        T_g * Dt(I_pz) ~ Ipcmd_in.u -  I_pz
        I_lvpl ~ LVPLogic(V, Zerox, Brkpt, L_vpl1)
        L_vplsw_adjust ~ ifelse(V>Brkpt, 0, L_vplsw)
        # Apply LVPL limiting and rate limiting in one step
        I_p_target ~ ifelse(L_vplsw_adjust > 0.5, uplimit(I_pz, I_lvpl), I_pz)
        Dt(I_p) ~ limiter(I_p_target - I_p, -rrpwr, rrpwr)
        #outputs
        Iqout.u ~ I_q
        Ipout.u ~ I_p
    end
end


