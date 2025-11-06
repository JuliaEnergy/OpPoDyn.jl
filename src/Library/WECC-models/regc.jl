@mtkmodel regc_a begin
    @structural_parameters begin
        L_vplsw=false
    end
    @components begin
        # inputs
        Vt_in = RealInput(guess=1)
        Iqcmd_in = RealInput(guess=-0.0567)
        Ipcmd_in = RealInput(guess=0.015)
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
        lvpnt0, [description=""]
        lvpnt1, [description=""]
        L_vpl1, [description="LVPL gain breakpoint (pu current on mbase / pu voltage)"]
        rrpwr, [description="Active current up-ramp rate limit on voltage recovery (pu/s)"]
        V_0lim, [description=""]
        K_hv, [description=""]
        I_0lim, [description=""]
    end
    @variables begin
        I_qrsum(t), [guess=0, description=""]
        I_qrlim(t), [guess=0, description=""]
        I_qr(t), [guess=0.0567, description=""]
        ΔV(t), [guess=-0.2, description=""]
        I_hv(t), [guess=-0.14, description=""]
        I_hvlim(t), [guess=0, description=""]
        I_q(t), [guess=0.0567, description="I_q after inverter current regulator with rate limits"]
        ΔI_q(t), [guess=0.0567, description=""]
        ΔI_pr(t), [guess=0, description=""]
        I_pr(t), [guess=0.015, description=""]
        ΔI_prlim(t), [guess=0, description=""]
        I_pg(t), [guess=0.015, description=""]
        y(t), [guess=1, description=""]
        I_p(t), [guess=0.015, description="I_p after inverter current regulator with limits"]
        V(t), [guess=1, description="V_t after filter"]
        I_lvpl(t), [guess=1.22, description="current resulting from low voltage power logic"]
    end
    @equations begin
        #q-phase current
        I_qrsum ~ -Iqcmd_in.u - I_qr
        I_qrlim ~ limiter(I_qrsum, I_qrmin, I_qrmax)
        T_g * Dt(I_qr) ~ I_qrlim
        ΔV ~ Vt_in.u - V_0lim
        I_hv ~ K_hv * ΔV
        I_hvlim ~ lowlimit(I_hv, 0)
        ΔI_q ~ I_qr - I_hvlim
        I_q ~ lowlimit(ΔI_q, I_0lim)
        #p-phase current
        ΔI_pr ~ Ipcmd_in.u - I_pr
        ΔI_prlim ~ uplimit(ΔI_pr, rrpwr)
        T_g * Dt(I_pg) ~ ΔI_prlim
        I_pr ~ ifelse(L_vplsw, uplimit(I_pg, I_lvpl), I_pg) #L_vplsw * uplimit(I_pg, I_lvpl) + (1-L_vplsw) * I_pg
        T_fltr * Dt(V) ~ Vt_in.u - V
        I_lvpl ~ LVPLogic(V, Zerox, Brkpt, L_vpl1)
        y ~ ifelse(Vt_in.u <= lvpnt0, 0, ifelse(Vt_in.u >= lvpnt1, 1, (Vt_in.u-lvpnt0)/(lvpnt1-lvpnt0)))
        I_p ~ y * I_pr
        #outputs
        Iqout.u ~ I_q
        Ipout.u ~ I_p
    end
end


