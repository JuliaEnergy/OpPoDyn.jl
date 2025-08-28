@mtkmodel reec_b begin
    @structural_parameters begin
        #vf_input = true
    end
    @components begin
        #terminal=Terminal()
        # inputs
        #if vf_input
        #    vf_in = RealInput(guess=1) # field voltage input [pu]
        #end
        Vt_in = RealInput(guess=1)
        P_e = RealInput(guess=1) #Inverter active power (pu on mbase)
        P_faref = RealInput(guess=1) #Inverter initial power factor angle (from power flow solution)
        Qext_in = RealInput(guess=1)
        Pref_in = RealInput(guess=1)
        Q_gen = RealInput(guess=1) #woher? Power Flow solution?
        # outputs
        Iqcmd_out = RealOutput()
        Ipcmd_out = RealOutput()
    end
    @parameters begin
        V_dip, [description="Low voltage condition trigger voltage (pu)"]
        V_up, [description="High voltage condition trigger voltage (pu)"]
        T_rv, [description="Terminal bus voltage filter time constant (s)"]
        V_ref0, [description="Reference voltage for reactive current injection (pu)"]
        dbd1, [description="Overvoltage deadband for reactive current injection (pu)"]
        dbd2, [description="Undervoltage deadband for reactive current injection (pu)"]
        K_qv, [description="Reactive current injection gain (pu/pu)"]
        I_qhl, [description="Maximum reactive current injection (pu on mbase)"]
        I_qll, [description="Minimum reactive current injection (pu on mbase)"]
        T_p, [description="Active power filter time constant (s)"]
        PfFlag, [description="Constant Q (0) or PF (1) local control"]
        Q_min, [description="Minimum reactive power when Vflag = 1 (pu on mbase)"]
        Q_max, [description="Maximum reactive power when Vflag = 1 (pu on mbase)"]
        V_min, [description="Minimum voltage at inverter terminal bus (pu)"]
        V_max, [description="Maximum voltage at inverter terminal bus (pu)"]
        K_qp, [description="Local Q regulator proportional gain (pu/pu)"]
        K_qi, [description="Local Q regulator integral gain (pu/pu-s)"]
        Vflag, [description="Local Q (0) or voltage control (1)"]
        K_vp, [description="Local voltage regulator proportional gain (pu/pu)"]
        K_vi, [description="Local voltage regulator integral gain (pu/pu-s)"]
        QFlag, [description="Bypass (0) or engage (1) inner voltage regulator loop"]
        I_max, [description="Maximum apparent current (pu on mbase)"]
        PqFlag, [description="Priority to reactive current (0) or active current (1)"]
        T_iq, [description="Reactive current regulator lag time constant (s)"]
        T_pord, [description="Inverter power order lag time constant (s)"]
        P_min, [description="Minimum active power (pu on mbase)"]
        P_max, [description="Maximum active power (pu on mbase)"]
        dP_min, [description="Active power down-ramp limit (pu/s on mbase)"]
        dP_max, [description="Active power up-ramp limit (pu/s on mbase)"]
    end
    @variables begin
        Voltage_dip(t), [description="freeze states if Voltagedip=1"]
        V_tfilt(t), [description="Voltage after filter"]
        ΔV_t(t), [description="Difference between filterd terminal voltage and reference voltage"]
        ΔV_tdbd(t), [description="Voltage after deadband"]
        I_qinj(t), [description="Limited current injection q-Phase from Voltage"]
        P_PF(t), [description="Inverter active power after filter"]
        Q_con(t), [description="Reactive Power after PfFlag"]
        Q_lim(t), [description="ReactivePower after limiter"]
        ΔQ(t), [description="Difference between Q_lim and Q_gen"]
        V_in(t), [description="Voltage after local Q regulator"]
        s_Q(t), [description="Frozen state in Q regulator"]
        V_lima(t), [description="Limited voltage after Q regulator"]
        V_con(t), [description="Voltage after Vflag"]
        V_limb(t), [description="Limited voltage V_con"]
        ΔV(t), [description="Difference between V_limb and V_tfilt"]
        I_in(t), [description="Current after local voltage regulator"]
        s_V(t), [description="Frozen state at local voltage regulator"]
        I_lim(t), [description="limited current after voltage regulator"]
        V_tfiltlim(t), [description="Voltage after filter with lower limit 0.01"]
        I_t(t), [description="Current from Q_con/V_tfiltlim"]
        I_qin(t), [description="Current after Reactive current regulator"]
        s_I(t), [description="Frozen state in reactive current regulator"]
        I_qcon(t), [description="Current after QFlag"]
        I_sum(t), [description="sum of I_qcon and I_qinj"]
        I_qcmd(t), [description="q-Phase output current"]
        P_refout(t), [description="Active power after inverter power order"]
        s_P(t), [description="frozen state after inverter power order"]
        P_lim(t), [description="Limited active power after inverter power order"]
        I_pref(t), [description="Current from P_lim/V_tfiltlim"]
        I_pcmd(t), [description="p-Phase output current"]
        I_qmin(t), [description="Minumum q-Phase current limit (pu)"]
        I_qmax(t), [description="Maximum q-Phase current limit (pu)"]
        I_pmax(t), [description="Maximum p-Phase current limit (pu)"]
        I_pmin(t), [description="Minumum p-Phase current limit (pu)"]
    end
    @equations begin
        Voltage_dip ~ ifelse(Vt_in.u<V_dip, 1, ifelse(Vt_in.u>V_up, 1, 0))
        T_rv * Dt(V_tfilt) ~ Vt_in.u - V_tfilt
        V_tfiltlim ~ ifelse(V_tfilt<0.01, 0.01, V_tfilt)
        #q-phase current
        ΔV_t ~ V_ref0 - V_tfilt
        ΔV_tdbd ~ deadband(ΔV_t, dbd1, dbd2)
        I_qinj ~ limiter(K_qp*ΔV_tdbd, I_qll, I_qhl)
        T_p * Dt(P_PF) ~ P_e.u - P_PF
        Q_con ~ PfFlag * P_PF * tan(P_faref.u) + (1-PfFlag) * Qext_in.u
        Q_lim ~ limiter(Q_con, Q_min, Q_max)
        ΔQ ~ Q_lim - Q_gen.u
        V_in ~ K_qp * ΔQ + s_Q
        Dt(s_Q) ~ (1-Voltage_dip) * K_qi * ΔQ
        V_lima ~ limiter(V_in, V_min, V_max)
        V_con ~ Vflag * V_lima + (1-Vflag) * Q_con
        V_limb ~ limiter(V_con, V_min, V_max)
        ΔV ~ V_limb - V_tfilt
        I_in ~ K_vp * ΔV + s_V
        Dt(s_V) ~ (1-Voltage_dip) * K_vi * ΔV
        I_lim ~ limiter(I_in, I_qmin, I_qmax)
        I_t ~ Q_con / V_tfiltlim
        T_iq * I_qin ~ s_I
        Dt(s_I) ~ (1-Voltage_dip) * (I_t - I_qin)
        I_qcon ~ QFlag * I_lim + (1-QFlag) * I_qin
        I_sum ~ I_qcon + I_qinj
        I_qcmd ~ limiter(I_sum, I_qmin, I_qmax)
        #p-phase current
        T_pord * P_refout ~ s_P
        Dt(s_P) ~ (1-Voltage_dip) * (Pref_in.u - P_refout)
        # P_lim with rate limiting - only define the derivative
        Dt(P_lim) ~ limiter(limiter(P_refout, P_min, P_max) - P_lim, dP_min, dP_max)
        I_pref ~ P_lim/V_tfiltlim
        I_pcmd ~ limiter(I_pref, I_pmin, I_pmax)
        #current limiter logic
        I_pmin ~ 0
        I_qmin ~ - I_qmax
        I_pmax ~ PqFlag * I_max + (1-PqFlag) * sqrt(I_max^2 - I_qcmd^2)
        I_qmax ~ PqFlag * sqrt(I_max^2 - I_pcmd^2) + (1-PqFlag) * I_max
        #outputs
        Iqcmd_out.u ~ I_qcmd
        Ipcmd_out.u ~ I_pcmd
    end
end



@mtkmodel reec_c begin
    @components begin
        Vt_in = RealInput(guess=1)
        P_e = RealInput(guess=1) #Inverter active power (pu on mbase)
        P_faref = RealInput(guess=1) #Inverter initial power factor angle (from power flow solution)
        Qext_in = RealInput(guess=1)
        Pref_in = RealInput(guess=1)
        P_gen = RealInput(guess=1)
        Q_gen = RealInput(guess=1) #woher? Power Flow solution?
        SOC = RealInput(guess=1)
        # outputs
        Iqcmd_out = RealOutput()
        Ipcmd_out = RealOutput()
    end
    @parameters begin
        V_dip, [description="voltage below which vltage-dip logic is initiated (pu)"]
        V_up, [description="voltage above which voltage dip/up logic is initiated(pu)"]
        T_rv, [description="voltage measurement transducer time constant"]
        V_ref0, [description="typically set to 1 pu (nominal voltage)"]
        dbd1, [description="lower voltage deadband error"]
        dbd2, [description="upper voltage deadband error"]
        K_qv, [description="Reactive current injection proportional gain (0 to disable)"]
        I_qh1, [description="Maximum reactive current injection (pu on mbase)"]
        I_ql1, [description="Minimum reactive current injection (pu on mbase)"]
        SOCini, [description="initial state of charge; must be between SOCmin and SOCmax"]
        SOCmax, [description="Maximum allowable state of charge in pu"]
        SOCmin, [description="Minimum allowable state of charge in pu"]
        T_p, [description="Active power filter time constant (s)"]
        PfFlag, [description="Constant Q (0) or PF (1) local control"]
        Q_min, [description="Minimum reactive power when Vflag = 1 (pu on mbase)"]
        Q_max, [description="Maximum reactive power when Vflag = 1 (pu on mbase)"]
        V_min, [description="Minimum voltage at inverter terminal bus (pu)"]
        V_max, [description="Maximum voltage at inverter terminal bus (pu)"]
        K_qp, [description="Local Q regulator proportional gain (pu/pu)"]
        K_qi, [description="Local Q regulator integral gain (pu/pu-s)"]
        Vflag, [description="Local Q (0) or voltage control (1)"]
        K_vp, [description="Local voltage regulator proportional gain (pu/pu)"]
        K_vi, [description="Local voltage regulator integral gain (pu/pu-s)"]
        QFlag, [description="Bypass (0) or engage (1) inner voltage regulator loop"]
        I_max, [description="Maximum apparent current (pu on mbase)"]
        PqFlag, [description="Priority to reactive current (0) or active current (1)"]
        T_iq, [description="Reactive current regulator lag time constant (s)"]
        T_pord, [description="Inverter power order lag time constant (s)"]
        P_min, [description="Minimum active power (pu on mbase)"]
        P_max, [description="Maximum active power (pu on mbase)"]
        dP_min, [description="Active power down-ramp limit (pu/s on mbase)"]
        dP_max, [description="Active power up-ramp limit (pu/s on mbase)"]
    end
    @variables begin
        Voltage_dip(t), [description="freeze states if Voltagedip=1"]
        V_tfilt(t), [description="Voltage after filter"]
        ΔV_t(t), [description="Difference between filterd terminal voltage and reference voltage"]
        ΔV_tdbd(t), [description="Voltage after deadband"]
        I_qinj(t), [description="Limited current injection q-Phase from Voltage"]
        P_PF(t), [description="Inverter active power after filter"]
        Q_con(t), [description="Reactive Power after PfFlag"]
        Q_lim(t), [description="ReactivePower after limiter"]
        ΔQ(t), [description="Difference between Q_lim and Q_gen"]
        V_in(t), [description="Voltage after local Q regulator"]
        s_Q(t), [description="Frozen state in Q regulator"]
        V_lima(t), [description="Limited voltage after Q regulator"]
        V_con(t), [description="Voltage after Vflag"]
        V_limb(t), [description="Limited voltage V_con"]
        ΔV(t), [description="Difference between V_limb and V_tfilt"]
        I_in(t), [description="Current after local voltage regulator"]
        s_V(t), [description="Frozen state at local voltage regulator"]
        I_lim(t), [description="limited current after voltage regulator"]
        V_tfiltlim(t), [description="Voltage after filter with lower limit 0.01"]
        I_t(t), [description="Current from Q_con/V_tfiltlim"]
        I_qin(t), [description="Current after Reactive current regulator"]
        s_I(t), [description="Frozen state in reactive current regulator"]
        I_qcon(t), [description="Current after QFlag"]
        I_sum(t), [description="sum of I_qcon and I_qinj"]
        I_qcmd(t), [description="q-Phase output current"]
        P_refout(t), [description="Active power after inverter power order"]
        s_P(t), [description="frozen state after inverter power order"]
        P_lim(t), [description="Limited active power after inverter power order"]
        I_pref(t), [description="Current from P_lim/V_tfiltlim"]
        I_pcmd(t), [description="p-Phase output current"]
        I_qmin(t), [description="Minumum q-Phase current limit (pu)"]
        I_qmax(t), [description="Maximum q-Phase current limit (pu)"]
        I_pmax(t), [description="Maximum p-Phase current limit (pu)"]
        I_pmin(t), [description="Minumum p-Phase current limit (pu)"]
    end
    @equations begin
        Voltage_dip ~ ifelse(Vt_in.u<V_dip, 1, ifelse(Vt_in.u>V_up, 1, 0))
        T_rv * Dt(V_tfilt) ~ Vt_in.u - V_tfilt
        V_tfiltlim ~ ifelse(V_tfilt<0.01, 0.01, V_tfilt)
        #q-phase current
        ΔV_t ~ V_ref0 - V_tfilt
        ΔV_tdbd ~ deadband(ΔV_t, dbd1, dbd2)
        I_qinj ~ limiter(K_qp*ΔV_tdbd, I_qll, I_qhl)
        T_p * Dt(P_PF) ~ P_e.u - P_PF
        Q_con ~ PfFlag * P_PF * tan(P_faref.u) + (1-PfFlag) * Qext_in.u
        Q_lim ~ limiter(Q_con, Q_min, Q_max)
        ΔQ ~ Q_lim - Q_gen.u
        V_in ~ K_qp * ΔQ + s_Q
        Dt(s_Q) ~ (1-Voltage_dip) * K_qi * ΔQ
        V_lima ~ limiter(V_in, V_min, V_max)
        V_con ~ Vflag * V_lima + (1-Vflag) * Q_con
        V_limb ~ limiter(V_con, V_min, V_max)
        ΔV ~ V_limb - V_tfilt
        I_in ~ K_vp * ΔV + s_V
        Dt(s_V) ~ (1-Voltage_dip) * K_vi * ΔV
        I_lim ~ limiter(I_in, I_qmin, I_qmax)
        I_t ~ Q_con / V_tfiltlim
        T_iq * I_qin ~ s_I
        Dt(s_I) ~ (1-Voltage_dip) * (I_t - I_qin)
        I_qcon ~ QFlag * I_lim + (1-QFlag) * I_qin
        I_sum ~ I_qcon + I_qinj
        I_qcmd ~ limiter(I_sum, I_qmin, I_qmax)
        #p-phase current
        T_pord * P_refout ~ s_P
        Dt(s_P) ~ (1-Voltage_dip) * (Pref_in.u - P_refout)
        # P_lim with rate limiting - only define the derivative
        Dt(P_lim) ~ limiter(limiter(P_refout, P_min, P_max) - P_lim, dP_min, dP_max)
        I_pref ~ P_lim/V_tfiltlim
        I_pcmd ~ limiter(I_pref, I_pmin, I_pmax)
        #current limiter logic
        I_pmin ~ 0
        I_qmin ~ - I_qmax
        I_pmax ~ PqFlag * I_max + (1-PqFlag) * sqrt(I_max^2 - I_qcmd^2)
        I_qmax ~ PqFlag * sqrt(I_max^2 - I_pcmd^2) + (1-PqFlag) * I_max
        #outputs
        Iqcmd_out.u ~ I_qcmd
        Ipcmd_out.u ~ I_pcmd
    end
end