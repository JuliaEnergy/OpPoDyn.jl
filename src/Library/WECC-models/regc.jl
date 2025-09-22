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
        lvpnt0, [description=""]
        lvpnt1, [description=""]
        L_vpl1, [description="LVPL gain breakpoint (pu current on mbase / pu voltage)"]
        L_vplsw, [description="Enable (1) or disable (0) low voltage power logic"]
        rrpwr, [description="Active current up-ramp rate limit on voltage recovery (pu/s)"]
        V_0lim, [description=""]
        K_hv, [description=""]
        I_0lim, [description=""]
    end
    @variables begin
        #I_qz(t), [description="I_q after inverter current regulator"]
        I_qrsum(t), [description=""]
        I_qrlim(t), [description=""]
        I_qr(t), [guess=1, description=""]
        ΔV(t), [description=""]
        I_hv(t), [description=""]
        I_hvlim(t), [description=""]
        I_q(t), [description="I_q after inverter current regulator with rate limits"]
        ΔI_q(t), [description=""]
        ΔI_pr(t), [description=""]
        I_pr(t), [description=""]
        ΔI_prlim(t), [description=""]
        I_pg(t), [description=""]
        y(t), [description=""]
        #I_pz(t), [description="I_p after Inverter current regulator"]
        I_p(t), [description="I_p after inverter current regulator with limits"]
        V(t), [description="V_t after filter"]
        I_lvpl(t), [description="current resulting from low voltage power logic"]
        #L_vplsw_adjust(t), [description="switch low voltage power logic off if V > breakpoint"]
    end
    @equations begin
        #q-phase current
        #T_g * Dt(I_qz) ~ -Iqcmd_in.u - I_qz
        #Dt(I_q) ~ ifelse(Q_gen0.u>0, uplimit(Dt(I_qz), I_qrmax), lowlimit(Dt(I_qz), I_qrmin))
        I_qrsum ~ -Iqcmd_in.u - I_qr
        I_qrlim ~ limiter(I_qrsum, I_qrmin, I_qrmax)
        T_g * Dt(I_qr) ~ I_qrlim
        ΔV ~ Vt_in.u - V_0lim
        I_hv ~ K_hv * ΔV
        I_hvlim ~ lowlimit(I_hv, 0)
        ΔI_q ~ I_qr - I_hvlim
        I_q ~ lowlimit(ΔI_q, I_0lim)
        #p-phase current
        #T_fltr * Dt(V) ~ Vt_in.u - V
        #T_g * Dt(I_pz) ~ Ipcmd_in.u -  I_pz
        #I_lvpl ~ LVPLogic(V, Zerox, Brkpt, L_vpl1)
        #L_vplsw_adjust ~ ifelse(V>Brkpt, 0, L_vplsw)
        #L_vplsw_adjust * I_p ~ L_vplsw_adjust * uplimit(I_pz, I_lvpl) #zwei Gleichungen für I_p. Wie ist das gemeint?
        #Dt(I_p) ~ uplimit(Dt(I_pz), rrpwr)
        ΔI_pr ~ Ipcmd_in.u - I_pr
        ΔI_prlim ~ uplimit(ΔI_pr, rrpwr)
        T_g * Dt(I_pg) ~ ΔI_prlim
        I_pr ~ L_vplsw * uplimit(I_pg, L_vpl1) + (1-L_vplsw) * I_pg
        T_fltr * Dt(V) ~ Vt_in.u - V
        I_lvpl ~ LVPLogic(V, Zerox, Brkpt, L_vpl1)
        y ~ ifelse(Vt_in.u > lvpnt0, ifelse(Vt_in.u < lvpnt1, (Vt_in.u-lvpnt0)/(lvpnt1-lvpnt0),1),0)
        I_p ~ y * I_pr
        #outputs
        Iqout.u ~ I_q
        Ipout.u ~ I_p
        #Iqout.u ~ I_hvlim #Iqcmd_in.u
        #Ipout.u ~ Ipcmd_in.u
    end
end


