@mtkmodel regc_a begin
    @structural_parameters begin
        L_vplsw=false
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
    @components begin
        # inputs
        Vt_in = RealInput(guess=1)
        Iqcmd_in = RealInput(guess=-0.056656797)
        Ipcmd_in = RealInput(guess=0.015)
        # outputs
        Iqout = RealOutput(guess=0.056656797)
        Ipout = RealOutput(guess=0.015000018)

        simpleLag = PowerDynamics.Library.SimpleLag(K=1, T=T_fltr, guess=1)
    end
    @variables begin
        I_qrsum(t), [guess=1.9317881e-14, description=""]
        I_qrlim(t), [guess=1.9317881e-14, description=""]
        I_qr(t), [guess=0.056656797, description=""]
        ΔV(t), [guess=-0.2, description=""]
        I_hv(t), [guess=-0.14, description=""]
        I_hvlim(t), [guess=0, description=""]
        I_q(t), [guess=0.056656797, description="I_q after inverter current regulator with rate limits"]
        ΔI_q(t), [guess=0.056656797, description=""]
        ΔI_pr(t), [guess=-7.359854e-12, description=""]
        I_pr(t), [guess=0.015, description=""]
        ΔI_prlim(t), [guess=-7.359854e-12, description=""]
        I_pg(t), [guess=0.015, description=""]
        y(t), [guess=1, description=""]
        I_p(t), [guess=0.015, description="I_p after inverter current regulator with limits"]
        V(t), [guess=1, description="V_t after filter"]
        I_lvpl(t), [guess=1.22, description="current resulting from low voltage power logic"]
    end
    @equations begin
        #q-phase current
        I_qrsum ~ -Iqcmd_in.u - I_qr
        I_qrlim ~ clamp(I_qrsum, I_qrmin, I_qrmax) #limiter(I_qrsum, I_qrmin, I_qrmax)
        T_g * Dt(I_qr) ~ I_qrlim
        ΔV ~ Vt_in.u - V_0lim
        I_hv ~ K_hv * ΔV
        I_hvlim ~ max(I_hv, 0) #lowlimit(I_hv, 0)
        ΔI_q ~ I_qr - I_hvlim
        I_q ~ max(ΔI_q, I_0lim) #lowlimit(ΔI_q, I_0lim)
        #p-phase current
        ΔI_pr ~ Ipcmd_in.u - I_pr
        ΔI_prlim ~ min(ΔI_pr, rrpwr) #uplimit(ΔI_pr, rrpwr)
        T_g * Dt(I_pg) ~ ΔI_prlim
        I_pr ~ ifelse(L_vplsw, min(I_pg, I_lvpl), I_pg)  #ifelse(L_vplsw, uplimit(I_pg, I_lvpl), I_pg) #L_vplsw * uplimit(I_pg, I_lvpl) + (1-L_vplsw) * I_pg

        #T_fltr * Dt(V) ~ Vt_in.u - V
        simpleLag.in ~ Vt_in.u
        V ~ simpleLag.out

        I_lvpl ~ LVPLogic(V, Zerox, Brkpt, L_vpl1)
        y ~ ifelse(Vt_in.u <= lvpnt0, 0, ifelse(Vt_in.u >= lvpnt1, 1, (Vt_in.u-lvpnt0)/(lvpnt1-lvpnt0)))
        I_p ~ y * I_pr
        #outputs
        Iqout.u ~ I_q
        Ipout.u ~ I_p
    end
end


@mtkmodel regc_a_pf begin
    @structural_parameters begin
        L_vplsw=false
    end
    @parameters begin
        I_qrmax, [description="Maximum rate-of-change of reactive current (pu/s)"]
        I_qrmin, [description="Minimum rate-of-change of reactive current (pu/s)"]
        T_gp, [description="Inverter current regulator lag time constant for active current (s)"]
        T_gq, [description="Inverter current regulator lag time constant for reactive current (s)"]
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
    @components begin
        # inputs
        #Vt_in = RealInput(guess=1)
        V_tfilt = RealInput(guess=1.001047) #TODO von außen eingeben! (aus reec_b) oder manuell hier berechnen?
        Iqcmd_in = RealInput(guess=-0.299956)
        Ipcmd_in = RealInput(guess=0.799884)
        # outputs
        Iqout = RealOutput(guess=0.299956)
        Ipout = RealOutput(guess=0.799884)

        simpleLag = PowerDynamics.Library.SimpleLag(K=1, T=T_fltr, guess=1.001047)
        SimpleLagLim = PowerDynamics.Library.SimpleLagLim(K=1, T=T_gq, outMin=I_qrmin, outMax=I_qrmax, guess=0.299956)
        SimpleLag_2uplims = PowerDynamics.Library.SimpleLag_2MaxLims(K=1, T=T_gp, doutMax=rrpwr, guess=0.799884, guessin=0.799884, guessx=0.799884)
    end
    @variables begin
        #I_qrsum(t), [guess=1.9317881e-14, description=""]
        #I_qrlim(t), [guess=1.9317881e-14, description=""]
        I_qr(t), [guess=0.299956, description=""]
        o2(t), [guess=0.30027, description=""]
        Hi_V_flag(t), [guess=0, description="If Vt<=V_0lim -> 0 ;  if Vt>V_0lim -> 1"]
        o3(t), [guess=0, description=""]
        V_tfiltlim(t), [guess=1.001047, description=""]
        ΔV(t), [guess=-0.198954, description=""]
        I_hv(t), [guess=-0.139268, description=""]
        I_hvlim(t), [guess=0, description=""]
        Q_gen(t), [guess=0.30027, description=""]
        I_q(t), [guess=0.299956, description="I_q after inverter current regulator with rate limits"]
        ΔI_q(t), [guess=0.30027, description=""]
        #ΔI_pr(t), [guess=-7.359854e-12, description=""]
        I_pr(t), [guess=0.799884, description=""]
        #ΔI_prlim(t), [guess=-7.359854e-12, description=""]
        #I_pg(t), [guess=0.015, description=""]
        Vt_scaled(t), [guess=1.001047, description=""]
        P_gen(t), [guess=0.800721, description=""]
        I_p(t), [guess=0.799884, description="I_p after inverter current regulator with limits"]
        V(t), [guess=1.001047, description="V_t after filter"]
        I_lvpl(t), [guess=999.99, description="current resulting from low voltage power logic"]
    end
    @equations begin
        #q-phase current
        SimpleLagLim.in ~ Iqcmd_in.u 
        I_qr ~ -SimpleLagLim.out #negative as in PF equation
        o2 ~ I_qr * V_tfilt.u
        Hi_V_flag ~ ifelse(V_tfilt.u>V_0lim, 1, 0)
        ΔV ~ V_tfilt.u - V_0lim
        I_hv ~ K_hv * ΔV
        I_hvlim ~ max(I_hv, 0)
        o3 ~ ifelse(Hi_V_flag==1, I_hvlim, 0)
        ΔI_q ~ o2 + o3
        Q_gen ~ max(ΔI_q, I_0lim) #TODO pf: Iolim1=Iolim*-1 und yo=min(yi,Iolim1)
        V_tfiltlim ~ max(V_tfilt.u, 0.01)
        I_q ~ Q_gen/V_tfiltlim

        #p-phase current
        #T_fltr * Dt(V) ~ Vt_in.u - V
        simpleLag.in ~ V_tfilt.u
        V ~ simpleLag.out

        I_lvpl ~ ifelse(L_vplsw, LVPLogic(V, Zerox, Brkpt, L_vpl1), 999.99) #no upper limit if switch of
        SimpleLag_2uplims.in ~ Ipcmd_in.u
        SimpleLag_2uplims.outMax ~ I_lvpl
        I_pr ~ SimpleLag_2uplims.out
        #P_gen ~ ifelse(Vt_in.u <= lvpnt0, 0, ifelse(Vt_in.u >= lvpnt1, 1, (Vt_in.u-lvpnt0)/(lvpnt1-lvpnt0)))
        Vt_scaled ~ V_tfilt.u * LVPLogic(V_tfilt.u, lvpnt0, lvpnt1, 1) #LVACM
        P_gen ~ Vt_scaled * I_pr
        I_p ~ P_gen/V_tfiltlim

        #outputs
        Iqout.u ~ I_q
        Ipout.u ~ I_p
    end
end


#Power Factory REGC_C WECC Generator-Converter Model
@mtkmodel regc_c_pf begin
    @structural_parameters begin
        #L_vplsw=false
    end
    @parameters begin
        I_qrmax, [description="Maximum rate-of-change of reactive current (pu/s)"]
        I_qrmin, [description="Minimum rate-of-change of reactive current (pu/s)"]
        #=
        T_gp, [description="Inverter current regulator lag time constant for active current (s)"]
        T_gq, [description="Inverter current regulator lag time constant for reactive current (s)"]
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
        =#
    end
    @components begin
        #=
        # inputs
        Vt_in = RealInput(guess=1) #TODO anscheinend gefilterters V_t schon! (PT1 Glied aus reec_b(_1))
        V_tfiltlim = RealInput(guess=1) #TODO von außen eingeben! (aus reec_b) oder manuell hier berechnen?
        Iqcmd_in = RealInput(guess=-0.056656797)
        Ipcmd_in = RealInput(guess=0.015)
        # outputs
        Iqout = RealOutput(guess=0.056656797)
        Ipout = RealOutput(guess=0.015000018)

        simpleLag = PowerDynamics.Library.SimpleLag(K=1, T=T_fltr, guess=1)
        SimpleLagLim = PowerDynamics.Library.SimpleLagLim(K=1, T=T_gq, outMin=I_qrmin, outMax=I_qrmax, guess=0.056656797)
        SimpleLag_2uplims = PowerDynamics.Library.SimpleLag_2MaxLims(K=1, T=T_gp, doutMax=rrpwr, guess=0.015)
        =#
        Iqcmd_in = RealInput(guess=-0.056656797)
        Ipcmd_in = RealInput(guess=0.015)
        Q_gen0 = RealInput(guess=0)
    end
    @variables begin
        #=
        I_qrsum(t), [guess=1.9317881e-14, description=""]
        I_qrlim(t), [guess=1.9317881e-14, description=""]
        I_qr(t), [guess=0.056656797, description=""]
        o2(t), [guess=0.056656797, description=""]
        Hi_V_flag(t), [guess=false, description="If Vt<=V_0lim -> 0 ;  if Vt>V_0lim -> 1"]
        o3(t), [guess=0, description=""]
        ΔV(t), [guess=-0.2, description=""]
        I_hv(t), [guess=-0.14, description=""]
        I_hvlim(t), [guess=0, description=""]
        Q_gen(t), [guess=0.056656797, description=""]
        I_q(t), [guess=0.056656797, description="I_q after inverter current regulator with rate limits"]
        ΔI_q(t), [guess=0.056656797, description=""]
        ΔI_pr(t), [guess=-7.359854e-12, description=""]
        I_pr(t), [guess=0.015, description=""]
        ΔI_prlim(t), [guess=-7.359854e-12, description=""]
        I_pg(t), [guess=0.015, description=""]
        Vt_scaled(t), [guess=1, description=""]
        P_gen(t), [guess=0.015, description=""]
        I_p(t), [guess=0.015, description="I_p after inverter current regulator with limits"]
        V(t), [guess=1, description="V_t after filter"]
        I_lvpl(t), [guess=1.22, description="current resulting from low voltage power logic"]
        =#
    end
    @equations begin
        I_qr ~ ifelse(Q_gen0.u<-0.000001, -clamp(I_qcmd.u, I_qrmin, 99999), ifelse(Q_gen0.u>0.000001, -clamp(I_qcmd.u, -99999, I_qrmax), -I_qcmd.u))


        #q-phase current
        SimpleLagLim.in ~ Iqcmd_in.u #TODO passt der Block -> -1 in PF? Block Definition hier anschauen!
        I_qr ~ SimpleLagLim.out
        o2 ~ I_qr * Vt_in.u
        Hi_V_flag ~ ifelse(Vt_in.u>V_0lim, true, false)
        ΔV ~ Vt_in.u - V_0lim
        I_hv ~ K_hv * ΔV
        I_hvlim ~ max(I_hv, 0)
        o3 ~ ifelse(Hi_V_flag, I_hvlim, 0)
        ΔI_q ~ o2 + o3
        Q_gen ~ max(ΔI_q, I_0lim) #lowlimit(ΔI_q, I_0lim) #TODO pf: Iolim1=Iolim*-1 und yo=min(yi,Iolim1)
        I_q ~ Q_gen/V_tfiltlim.u

        #p-phase current
        #T_fltr * Dt(V) ~ Vt_in.u - V
        simpleLag.in ~ Vt_in.u
        V ~ simpleLag.out

        I_lvpl ~ ifelse(L_vplsw, LVPLogic(V, Zerox, Brkpt, L_vpl1), 999) #no upper limit if switch of
        SimpleLag_2uplims.in ~ Ipcmd_in.u
        SimpleLag_2uplims.outMax ~ I_lvpl
        I_pr ~ SimpleLag_2uplims.out
        #P_gen ~ ifelse(Vt_in.u <= lvpnt0, 0, ifelse(Vt_in.u >= lvpnt1, 1, (Vt_in.u-lvpnt0)/(lvpnt1-lvpnt0)))
        Vt_scaled ~ Vt_in.u * LVPLogic(Vt_in.u, lvpnt0, lvpnt1, 1)
        P_gen ~ Vt_scaled * I_pr
        I_p ~ P_gen/V_tfiltlim.u

        #outputs
        Iqout.u ~ I_q
        Ipout.u ~ I_p
    end
end

