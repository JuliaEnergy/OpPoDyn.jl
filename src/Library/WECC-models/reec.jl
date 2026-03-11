@mtkmodel reec_b begin
    @structural_parameters begin
        PfFlag = false
        Vflag = false
        QFlag = false
        PqFlag = false
    end
    @parameters begin
        V_dip, [description="Low voltage condition trigger voltage (pu)"]
        V_up, [description="High voltage condition trigger voltage (pu)"]
        T_rv, [description="Terminal bus voltage filter time constant (s)"]
        V_ref0, [description="Reference voltage for reactive current injection (pu)"]
        dbd1, [description="Overvoltage deadband for reactive current injection (pu)"]
        dbd2, [description="Undervoltage deadband for reactive current injection (pu)"]
        K_qv, [description="Reactive current injection gain (pu/pu)"]
        I_qh1, [description="Maximum reactive current injection (pu on mbase)"]
        I_ql1, [description="Minimum reactive current injection (pu on mbase)"]
        T_p, [description="Active power filter time constant (s)"]
        Q_min, [description="Minimum reactive power when Vflag = 1 (pu on mbase)"]
        Q_max, [description="Maximum reactive power when Vflag = 1 (pu on mbase)"]
        V_min, [description="Minimum voltage at inverter terminal bus (pu)"]
        V_max, [description="Maximum voltage at inverter terminal bus (pu)"]
        K_qp, [description="Local Q regulator proportional gain (pu/pu)"]
        K_qi, [description="Local Q regulator integral gain (pu/pu-s)"]
        K_vp, [description="Local voltage regulator proportional gain (pu/pu)"]
        K_vi, [description="Local voltage regulator integral gain (pu/pu-s)"]
        I_max, [description="Maximum apparent current (pu on mbase)"]
        T_iq, [description="Reactive current regulator lag time constant (s)"]
        T_pord, [description="Inverter power order lag time constant (s)"]
        P_min, [description="Minimum active power (pu on mbase)"]
        P_max, [description="Maximum active power (pu on mbase)"]
        dP_min, [description="Active power down-ramp limit (pu/s on mbase)"]
        dP_max, [description="Active power up-ramp limit (pu/s on mbase)"]
    end
    @components begin
        Vt_in = RealInput(guess=1)
        P_e = RealInput(guess=0.015) #Inverter active power (pu on mbase)
        P_faref = RealInput(guess=-1.31199) #Inverter initial power factor angle (from power flow solution)
        Qext_in = RealInput(guess=-0.056656797)
        Pref_in = RealInput(guess=0.015)
        Q_gen = RealInput(guess=-0.056656801)
        # outputs
        Iqcmd_out = RealOutput(guess=-0.056656797)
        Ipcmd_out = RealOutput(guess=0.015)

        #building blocks
        simpleLag = PowerDynamics.Library.SimpleLag(K=1, T=T_rv, guess=1)
        simpleLag1 = PowerDynamics.Library.SimpleLag(K=1, T=T_p, guess=0.015)
        deadband = PowerDynamics.Library.DeadZone(uMax=dbd2, uMin=dbd1)
    end
    @variables begin
        Voltage_dip(t), [guess=0, description="freeze states if Voltagedip=1"]
        V_tfilt(t), [guess=1, description="Voltage after filter"]
        V_tfiltlim(t), [guess=1, description="Voltage after filter with lower limit 0.01"]
        ΔV_t(t), [guess=6.6613381e-14, description="Difference between filterd terminal voltage and reference voltage"]
        ΔV_tdbd(t), [guess=0, description="Voltage after deadband"]
        I_qinj(t), [guess=0, description="Limited current injection q-Phase from Voltage"]
        P_PF(t), [guess=0.015, description="Inverter active power after filter"]
        Q_con(t), [guess=-0.056656797, description="Reactive Power after PfFlag"]
        Q_lim(t), [guess=-0.056656797, description="ReactivePower after limiter"]
        ΔQ(t), [guess=4.8269116e-9, description="Difference between Q_lim and Q_gen"]
        s_Q(t), [guess=4.8269116e-9, description="Frozen state in Q regulator"]
        s_Qint(t), [guess=1, description=""]
        V_in(t), [guess=1, description="Voltage after local Q regulator"]
        V_lima(t), [guess=1, description="Limited voltage after Q regulator"]
        V_con(t), [guess=1, description="Voltage after Vflag"]
        V_limb(t), [guess=1, description="Limited voltage V_con"]
        ΔV(t), [guess=6.6613381e-14, description="Difference between V_limb and V_tfilt"]
        s_V(t), [guess=6.6613381e-14, description="Frozen state at local voltage regulator"]
        s_Vint(t), [guess=-0.0567, description=""]
        I_in(t), [guess=-0.0567, description="Current after local voltage regulator"]
        I_lim(t), [guess=-0.0567, description="limited current after voltage regulator"]
        I_t(t), [guess=-0.056656797, description="Current from Q_con/V_tfiltlim"]
        ΔI(t), [guess=-1.9310942e-14, description=""]
        I_qin(t), [guess=-0.056656797, description="Current after Reactive current regulator"]
        I_qcon(t), [guess=-0.056656797, description="Current after QFlag"]
        I_sum(t), [guess=-0.056656797, description="sum of I_qcon and I_qinj"]
        I_qcmd(t), [guess=-0.056656797, description="q-Phase output current"]
        P_refout(t), [guess=0.015, description="Active power after inverter power order"]
        P_lim(t), [guess=0.015, description="Limited active power after inverter power order"]
        ΔP(t), [guess=-7.3598558e-12, description="Active power difference between P_ref and P_refout"]
        ΔP_lim(t), [guess=-7.3598558e-12, description="Ramp-limited active power difference"]
        I_pref(t), [guess=0.015, description="Current from P_lim/V_tfiltlim"]
        I_pcmd(t), [guess=0.015, description="p-Phase output current"]
        I_qmin(t), [guess=-1.82, description="Minumum q-Phase current limit (pu)"]
        I_qmax(t), [guess=1.82, description="Maximum q-Phase current limit (pu)"]
        I_pmax(t), [guess=1.8191179, description="Maximum p-Phase current limit (pu)"]
        I_pmin(t), [guess=0, description="Minumum p-Phase current limit (pu)"]
    end
    @equations begin
        Voltage_dip ~ ifelse(Vt_in.u<V_dip, 1, ifelse(Vt_in.u>V_up, 1, 0)) #0

        #T_rv * Dt(V_tfilt) ~ Vt_in.u - V_tfilt
        simpleLag.in ~ Vt_in.u
        V_tfilt ~ simpleLag.out

        V_tfiltlim ~ max(V_tfilt, 0.01) #lowlimit(V_tfilt, 0.01)
        #q-phase current
        ΔV_t ~ V_ref0 - V_tfilt

        #ΔV_tdbd ~ deadband(ΔV_t, dbd1, dbd2)
        deadband.in ~ ΔV_t
        ΔV_tdbd ~ deadband.out

        I_qinj ~ clamp(K_qv*ΔV_tdbd, I_ql1, I_qh1) #limiter(K_qv*ΔV_tdbd, I_ql1, I_qh1)

        #T_p * Dt(P_PF) ~ P_e.u - P_PF
        simpleLag1.in ~ P_e.u
        P_PF ~ simpleLag1.out

        Q_con ~ ifelse(PfFlag, P_PF * tan(P_faref.u), Qext_in.u) #PfFlag * P_PF * tan(P_faref.u) + (1-PfFlag) * Qext_in.u
        Q_lim ~ clamp(Q_con, Q_min, Q_max) #limiter(Q_con, Q_min, Q_max)
        ΔQ ~ Q_lim - Q_gen.u
        s_Q ~ (1-Voltage_dip) * ΔQ
        Dt(s_Qint) ~ K_qi * s_Q
        V_in ~ K_qp * s_Q + s_Qint
        V_lima ~ clamp(V_in, V_min, V_max) #limiter(V_in, V_min, V_max)
        V_con ~ ifelse(Vflag, V_lima, V_ref0) #Vflag * V_lima + (1-Vflag) * V_ref0
        V_limb ~ clamp(V_con, V_min, V_max) #limiter(V_con, V_min, V_max)
        ΔV ~ V_limb - V_tfilt
        s_V ~ (1-Voltage_dip) * ΔV
        Dt(s_Vint) ~ K_vi * s_V
        I_in ~ K_vp * s_V + s_Vint
        I_lim ~ clamp(I_in, I_qmin, I_qmax) #limiter(I_in, I_qmin, I_qmax)
        I_t ~ Q_con / V_tfiltlim
        ΔI ~ I_t - I_qin
        T_iq * Dt(I_qin) ~ (1-Voltage_dip) * ΔI
        I_qcon ~ ifelse(QFlag, I_lim, I_qin) #QFlag * I_lim + (1-QFlag) * I_qin
        I_sum ~ I_qcon + I_qinj
        I_qcmd ~ clamp(I_sum, I_qmin, I_qmax) #limiter(I_sum, I_qmin, I_qmax)
        #p-phase current
        ΔP ~ Pref_in.u - P_refout
        ΔP_lim ~ clamp(ΔP, dP_min, dP_max) #limiter(ΔP, dP_min, dP_max)
        T_pord * Dt(P_refout) ~ (1-Voltage_dip) * ΔP_lim
        P_lim ~ clamp(P_refout, P_min, P_max) #limiter(P_refout, P_min, P_max)
        I_pref ~ P_lim/V_tfiltlim
        I_pcmd ~ clamp(I_pref, I_pmin, I_pmax) #limiter(I_pref, I_pmin, I_pmax)
        #current limiter logic
        I_pmin ~ 0
        I_qmin ~ - I_qmax
        I_pmax ~ ifelse(PqFlag, I_max, sqrt(I_max^2 - I_qcmd^2)) #PqFlag * I_max + (1-PqFlag) * sqrt(I_max^2 - I_qcmd^2)
        I_qmax ~ ifelse(PqFlag, sqrt(I_max^2 - I_pcmd^2), I_max) #PqFlag * sqrt(I_max^2 - I_pcmd^2) + (1-PqFlag) * I_max
        #outputs
        Iqcmd_out.u ~ I_qcmd
        Ipcmd_out.u ~ I_pcmd
    end
end



@mtkmodel reec_c begin
    @structural_parameters begin
        PfFlag = false
        Vflag = false
        QFlag = false
        PqFlag = false
    end
    @parameters begin
        V_dip, [description="Low voltage condition trigger voltage (pu)"]
        V_up, [description="High voltage condition trigger voltage (pu)"]
        T_rv, [description="Terminal bus voltage filter time constant (s)"]
        V_ref0, [description="Reference voltage for reactive current injection (pu)"]
        dbd1, [description="Overvoltage deadband for reactive current injection (pu)"]
        dbd2, [description="Undervoltage deadband for reactive current injection (pu)"]
        K_qv, [description="Reactive current injection gain (pu/pu)"]
        I_qh1, [description="Maximum reactive current injection (pu on mbase)"]
        I_ql1, [description="Minimum reactive current injection (pu on mbase)"]
        T_p, [description="Active power filter time constant (s)"]
        Q_min, [description="Minimum reactive power when Vflag = 1 (pu on mbase)"]
        Q_max, [description="Maximum reactive power when Vflag = 1 (pu on mbase)"]
        V_min, [description="Minimum voltage at inverter terminal bus (pu)"]
        V_max, [description="Maximum voltage at inverter terminal bus (pu)"]
        K_qp, [description="Local Q regulator proportional gain (pu/pu)"]
        K_qi, [description="Local Q regulator integral gain (pu/pu-s)"]
        K_vp, [description="Local voltage regulator proportional gain (pu/pu)"]
        K_vi, [description="Local voltage regulator integral gain (pu/pu-s)"]
        I_max, [description="Maximum apparent current (pu on mbase)"]
        T_iq, [description="Reactive current regulator lag time constant (s)"]
        T_pord, [description="Inverter power order lag time constant (s)"]
        P_min, [description="Minimum active power (pu on mbase)"]
        P_max, [description="Maximum active power (pu on mbase)"]
        dP_min, [description="Active power down-ramp limit (pu/s on mbase)"]
        dP_max, [description="Active power up-ramp limit (pu/s on mbase)"]
        soc_ini, [description="initial state of charge"]
        T_char, [description="Battery discharge time"]
        SOCmin, [description="Minimum allowable state of charge"]
        SOCmax, [description="Maximum allowable state of charge"]
        Vq1=0.0, [description="q-VDL Table"]
        Vq2=0.2, [description="q-VDL Table"]
        Vq3=0.5, [description="q-VDL Table"]
        Vq4=1, [description="q-VDL Table"]
        Iq1=0.75, [description="q-VDL Table"]
        Iq2=0.75, [description="q-VDL Table"]
        Iq3=0.75, [description="q-VDL Table"]
        Iq4=0.75, [description="q-VDL Table"]
        Vp1=0.2, [description="p-VDL Table"]
        Vp2=0.5, [description="p-VDL Table"]
        Vp3=0.75, [description="p-VDL Table"]
        Vp4=1, [description="p-VDL Table"]
        Ip1=1.11, [description="p-VDL Table"]
        Ip2=1.11, [description="p-VDL Table"]
        Ip3=1.11, [description="p-VDL Table"]
        Ip4=1.11, [description="p-VDL Table"]
    end
    @components begin
        Vt_in = RealInput(guess=1)
        P_e = RealInput(guess=0.015) #Inverter active power (pu on mbase)
        P_faref = RealInput(guess=-1.31199) #Inverter initial power factor angle (from power flow solution)
        Qext_in = RealInput(guess=-0.056656797)
        Pref_in = RealInput(guess=0.015)
        Q_gen = RealInput(guess=-0.056656801)
        P_aux = RealInput(guess=0)
        PELEC = RealInput(guess=0.015)
        # outputs
        Iqcmd_out = RealOutput(guess=-0.056656797)
        Ipcmd_out = RealOutput(guess=0.015)

        #building blocks
        simpleLag = PowerDynamics.Library.SimpleLag(K=1, T=T_rv, guess=1)
        simpleLag1 = PowerDynamics.Library.SimpleLag(K=1, T=T_p, guess=0.015)
        deadband = PowerDynamics.Library.DeadZone(uMax=dbd2, uMin=dbd1)
    end
    @variables begin
        Voltage_dip(t), [guess=0, description="freeze states if Voltagedip=1"]
        V_tfilt(t), [guess=1, description="Voltage after filter"]
        V_tfiltlim(t), [guess=1, description="Voltage after filter with lower limit 0.01"]
        ΔV_t(t), [guess=6.6613381e-14, description="Difference between filterd terminal voltage and reference voltage"]
        ΔV_tdbd(t), [guess=0, description="Voltage after deadband"]
        I_qinj(t), [guess=0, description="Limited current injection q-Phase from Voltage"]
        P_PF(t), [guess=0.015, description="Inverter active power after filter"]
        Q_con(t), [guess=-0.056656797, description="Reactive Power after PfFlag"]
        Q_lim(t), [guess=-0.056656797, description="ReactivePower after limiter"]
        ΔQ(t), [guess=4.8269116e-9, description="Difference between Q_lim and Q_gen"]
        s_Q(t), [guess=4.8269116e-9, description="Frozen state in Q regulator"]
        s_Qint(t), [guess=1, description=""]
        V_in(t), [guess=1, description="Voltage after local Q regulator"]
        V_lima(t), [guess=1, description="Limited voltage after Q regulator"]
        V_con(t), [guess=1, description="Voltage after Vflag"]
        V_limb(t), [guess=1, description="Limited voltage V_con"]
        ΔV(t), [guess=6.6613381e-14, description="Difference between V_limb and V_tfilt"]
        s_V(t), [guess=6.6613381e-14, description="Frozen state at local voltage regulator"]
        s_Vint(t), [guess=-0.0567, description=""]
        I_in(t), [guess=-0.0567, description="Current after local voltage regulator"]
        I_lim(t), [guess=-0.0567, description="limited current after voltage regulator"]
        I_t(t), [guess=-0.056656797, description="Current from Q_con/V_tfiltlim"]
        ΔI(t), [guess=0, description=""]
        I_qin(t), [guess=-0.056656797, description="Current after Reactive current regulator"]
        I_qcon(t), [guess=-0.056656797, description="Current after QFlag"]
        I_sum(t), [guess=-0.056656797, description="sum of I_qcon and I_qinj"]
        I_qcmd(t), [guess=-0.056656797, description="q-Phase output current"]
        P_refout(t), [guess=0.015, description="Active power after inverter power order"]
        P_lim(t), [guess=0.015, description="Limited active power after inverter power order"]
        ΔP(t), [guess=0, description="Active power difference between P_ref and P_refout"]
        ΔP_lim(t), [guess=0, description="Ramp-limited active power difference"]
        I_pref(t), [guess=0.015, description="Current from P_lim/V_tfiltlim"]
        ΔI_p(t), [guess=0.015, description=""]
        I_pcmd(t), [guess=0.015, description="p-Phase output current"]
        I_qmin(t), [guess=-0.75, description="Minumum q-Phase current limit (pu)"]
        I_qmax(t), [guess=0.75, description="Maximum q-Phase current limit (pu)"]
        I_pmax(t), [guess=1.108553114, description="Maximum p-Phase current limit (pu)"]
        I_pmin(t), [guess=-1.108553114, description="Minumum p-Phase current limit (pu)"]
        I_pmin_soc(t), [guess=-1.108553114, description=""]
        I_pmax_soc(t), [guess=1.108553114, description=""]
        soc_Imin(t), [guess=1, description=""]
        soc_Imax(t), [guess=1, description=""]
        P_stor(t), [guess=0.015, description=""]
        soc(t), [guess=0.485, description=""]
        soc_lim(t), [guess=0.485, description=""]
        VDL1_out(t), [guess=0.75, description=""]
        VDL2_out(t), [guess=1.11, description=""]
    end
    @equations begin
        Voltage_dip ~ ifelse(Vt_in.u<V_dip, 1, ifelse(Vt_in.u>V_up, 1, 0))

        #T_rv * Dt(V_tfilt) ~ Vt_in.u - V_tfilt
        simpleLag.in ~ Vt_in.u
        V_tfilt ~ simpleLag.out

        V_tfiltlim ~ max(V_tfilt, 0.01) #lowlimit(V_tfilt, 0.01)
        #q-phase current
        ΔV_t ~ V_ref0 - V_tfilt

        #ΔV_tdbd ~ deadband(ΔV_t, dbd1, dbd2)
        deadband.in ~ ΔV_t
        ΔV_tdbd ~ deadband.out

        I_qinj ~ clamp(K_qv*ΔV_tdbd, I_ql1, I_qh1) #limiter(K_qv*ΔV_tdbd, I_ql1, I_qh1)

        #T_p * Dt(P_PF) ~ P_e.u - P_PF
        simpleLag1.in ~ P_e.u
        P_PF ~ simpleLag1.out

        Q_con ~ ifelse(PfFlag, P_PF * tan(P_faref.u), Qext_in.u) #PfFlag * P_PF * tan(P_faref.u) + (1-PfFlag) * Qext_in.u
        Q_lim ~ clamp(Q_con, Q_min, Q_max) #limiter(Q_con, Q_min, Q_max)
        ΔQ ~ Q_lim - Q_gen.u
        s_Q ~ (1-Voltage_dip) * ΔQ
        Dt(s_Qint) ~ K_qi * s_Q
        V_in ~ K_qp * s_Q + s_Qint
        V_lima ~ clamp(V_in, V_min, V_max) #limiter(V_in, V_min, V_max)
        V_con ~ ifelse(Vflag, V_lima, V_ref0) #Vflag * V_lima + (1-Vflag) * V_ref0
        V_limb ~ clamp(V_con, V_min, V_max) #limiter(V_con, V_min, V_max)
        ΔV ~ V_limb - V_tfilt
        s_V ~ (1-Voltage_dip) * ΔV
        Dt(s_Vint) ~ K_vi * s_V
        I_in ~ K_vp * s_V + s_Vint
        I_lim ~ clamp(I_in, I_qmin, I_qmax) #limiter(I_in, I_qmin, I_qmax)
        I_t ~ Q_con / V_tfiltlim
        ΔI ~ I_t - I_qin
        T_iq * Dt(I_qin) ~ (1-Voltage_dip) * ΔI
        I_qcon ~ ifelse(QFlag, I_lim, I_qin) #QFlag * I_lim + (1-QFlag) * I_qin
        I_sum ~ I_qcon + I_qinj
        I_qcmd ~ clamp(I_sum, I_qmin, I_qmax) #limiter(I_sum, I_qmin, I_qmax)
        #p-phase current
        ΔP ~ Pref_in.u - P_refout
        ΔP_lim ~ clamp(ΔP, dP_min, dP_max) #limiter(ΔP, dP_min, dP_max)
        T_pord * Dt(P_refout) ~ (1-Voltage_dip) * ΔP_lim
        P_lim ~ clamp(P_refout, P_min, P_max) #limiter(P_refout, P_min, P_max)
        I_pref ~ P_lim/V_tfiltlim
        ΔI_p ~ P_aux.u + I_pref
        I_pmin_soc ~ I_pmin * soc_Imin
        I_pmax_soc ~ I_pmax * soc_Imax
        I_pcmd ~ clamp(ΔI_p, I_pmin_soc, I_pmax_soc) #limiter(ΔI_p, I_pmin_soc, I_pmax_soc)
        #soc logic
        T_char * Dt(P_stor) ~ PELEC.u
        soc ~ soc_ini - P_stor
        soc_lim ~ clamp(soc, SOCmin, SOCmax) #limiter(soc, SOCmin, SOCmax)
        soc_Imax ~ ifelse(soc_lim<=SOCmin, 0, 1)
        soc_Imin ~ ifelse(soc_lim>=SOCmax, 0, 1)
        #VDL tables
        VDL1_out ~ VDL(V_tfilt, Vq1, Vq2, Vq3, Vq4, Iq1, Iq2, Iq3, Iq4)
        VDL2_out ~ VDL(V_tfilt, Vp1, Vp2, Vp3, Vp4, Ip1, Ip2, Ip3, Ip4)
        #current limiter logic
        I_pmin ~ -I_pmax
        I_qmin ~ -I_qmax
        I_pmax ~ ifelse(PqFlag, min(VDL2_out, I_max), min(VDL2_out, sqrt(I_max^2 - I_qcmd^2))) #PqFlag * min(VDL2_out, I_max) + (1-PqFlag) * min(VDL2_out, sqrt(I_max^2 - I_qcmd^2))
        I_qmax ~ ifelse(PqFlag, min(VDL1_out, sqrt(I_max^2 - I_pcmd^2)), min(VDL1_out, I_max)) #PqFlag * min(VDL1_out, sqrt(I_max^2 - I_pcmd^2)) + (1-PqFlag) * min(VDL1_out, I_max)
        #outputs
        Iqcmd_out.u ~ I_qcmd
        Ipcmd_out.u ~ I_pcmd
    end
end


@mtkmodel reec_a begin
    @structural_parameters begin
        PfFlag = false
        Vflag = false
        QFlag = false
        PqFlag = false
    end
    @parameters begin
        V_0, [description=""]
        V_dip, [description="Low voltage condition trigger voltage (pu)"]
        V_up, [description="High voltage condition trigger voltage (pu)"]
        T_rv, [description="Terminal bus voltage filter time constant (s)"]
        V_ref0, [description="Reference voltage for reactive current injection (pu); has to be !=0"]
        dbd1, [description="Overvoltage deadband for reactive current injection (pu)"]
        dbd2, [description="Undervoltage deadband for reactive current injection (pu)"]
        K_qv, [description="Reactive current injection gain (pu/pu)"]
        I_qh1, [description="Maximum reactive current injection (pu on mbase)"]
        I_ql1, [description="Minimum reactive current injection (pu on mbase)"]
        T_p, [description="Active power filter time constant (s)"]
        Q_min, [description="Minimum reactive power when Vflag = 1 (pu on mbase)"]
        Q_max, [description="Maximum reactive power when Vflag = 1 (pu on mbase)"]
        V_min, [description="Minimum voltage at inverter terminal bus (pu)"]
        V_max, [description="Maximum voltage at inverter terminal bus (pu)"]
        K_qp, [description="Local Q regulator proportional gain (pu/pu)"]
        K_qi, [description="Local Q regulator integral gain (pu/pu-s)"]
        V_bias, [description=""]
        K_vp, [description="Local voltage regulator proportional gain (pu/pu)"]
        K_vi, [description="Local voltage regulator integral gain (pu/pu-s)"]
        I_max, [description="Maximum apparent current (pu on mbase)"]
        T_iq, [description="Reactive current regulator lag time constant (s)"]
        T_pord, [description="Inverter power order lag time constant (s)"]
        P_min, [description="Minimum active power (pu on mbase)"]
        P_max, [description="Maximum active power (pu on mbase)"]
        dP_min, [description="Active power down-ramp limit (pu/s on mbase)"]
        dP_max, [description="Active power up-ramp limit (pu/s on mbase)"]
        Vq1=0.1, [description="q-VDL Table 1"]
        Vq2=0.4, [description="q-VDL Table 1"]
        Vq3=0.6, [description="q-VDL Table 1"]
        Vq4=0.9, [description="q-VDL Table 1"]
        Iq1=0.01, [description="q-VDL Table 1"]
        Iq2=0.5, [description="q-VDL Table 1"]
        Iq3=0.7, [description="q-VDL Table 1"]
        Iq4=1.0, [description="q-VDL Table 1"]
        Vp1=0.1, [description="p-VDL Table 2"]
        Vp2=0.5, [description="p-VDL Table 2"]
        Vp3=0.9, [description="p-VDL Table 2"]
        Vp4=1, [description="p-VDL Table 2"]
        Ip1=0.4, [description="p-VDL Table 2"]
        Ip2=0.7, [description="p-VDL Table 2"]
        Ip3=1.2, [description="p-VDL Table 2"]
        Ip4=1.2, [description="p-VDL Table 2"]
    end
        @components begin
        Vt_in = RealInput(guess=1)
        P_e = RealInput(guess=0.015) #Inverter active power (pu on mbase)
        P_faref = RealInput(guess=-1.31199) #Inverter initial power factor angle (from power flow solution)
        Qext_in = RealInput(guess=-0.056656797)
        Pref_in = RealInput(guess=0.015)
        Q_gen = RealInput(guess=-0.056656801)
        Wg = RealInput(guess=1) #Rotational speed generator
        # outputs
        Iqcmd_out = RealOutput(guess=-0.056656797)
        Ipcmd_out = RealOutput(guess=0.015)

        #building blocks
        simpleLag = PowerDynamics.Library.SimpleLag(K=1, T=T_rv, guess=1)
        simpleLag1 = PowerDynamics.Library.SimpleLag(K=1, T=T_p, guess=0.015)
        deadband = PowerDynamics.Library.DeadZone(uMax=dbd2, uMin=dbd1)
    end
    @variables begin
        Voltage_dip(t), [guess=0, description="freeze states if Voltagedip=1"]
        V_tfilt(t), [guess=1, description="Voltage after filter"]
        V_tfiltlim(t), [guess=1, description="Voltage after filter with lower limit 0.01"]
        ΔV_t(t), [guess=0, description="Difference between filterd terminal voltage and reference voltage"]
        ΔV_tdbd(t), [guess=0, description="Voltage after deadband"]
        I_qinj(t), [guess=0, description="Limited current injection q-Phase from Voltage"]
        P_PF(t), [guess=0.015, description="Inverter active power after filter"]
        Q_con(t), [guess=-0.056656797, description="Reactive Power after PfFlag"]
        Q_lim(t), [guess=-0.056656797, description="ReactivePower after limiter"]
        ΔQ(t), [guess=4.8269116e-9, description="Difference between Q_lim and Q_gen"]
        s_Q(t), [guess=1, description="Frozen state in Q regulator"]
        V_in(t), [guess=1, description="Voltage after local Q regulator"]
        V_lima(t), [guess=1, description="Limited voltage after Q regulator"]
        V_mod(t), [guess=0, description=""]
        V_con(t), [guess=1, description="Voltage after Vflag"]
        V_limb(t), [guess=0.9, description="Limited voltage V_con"]
        ΔV(t), [guess=-0.1, description="Difference between V_limb and V_tfilt"]
        s_V(t), [guess=-0.056658, description="Frozen state at local voltage regulator"]
        I_in(t), [guess=-0.216658, description="Current after local voltage regulator"]
        I_lim(t), [guess=-0.216658, description="limited current after voltage regulator"]
        I_t(t), [guess=-0.056656797, description="Current from Q_con/V_tfiltlim"]
        ΔI(t), [guess=0, description=""]
        I_qin(t), [guess=-0.056656797, description="Current after Reactive current regulator"]
        I_qcon(t), [guess=-0.056656797, description="Current after QFlag"]
        I_sum(t), [guess=-0.056656797, description="sum of I_qcon and I_qinj"]
        I_qcmd(t), [guess=-0.056656797, description="q-Phase output current"]
        P_in(t), [guess=0.015, description=""]
        P_refout(t), [guess=0.015, description="Active power after inverter power order"]
        P_lim(t), [guess=0.015, description="Limited active power after inverter power order"]
        ΔP(t), [guess=0, description="Active power difference between P_ref and P_refout"]
        ΔP_lim(t), [guess=0, description="Ramp-limited active power difference"]
        I_pref(t), [guess=0.015, description="Current from P_lim/V_tfiltlim"]
        I_pcmd(t), [guess=0.015, description="p-Phase output current"]
        I_qmin(t), [guess=-1.1, description="Minumum q-Phase current limit (pu)"]
        I_qmax(t), [guess=1.1, description="Maximum q-Phase current limit (pu)"]
        I_pmax(t), [guess=0.9259688073, description="Maximum p-Phase current limit (pu)"]
        I_pmin(t), [guess=0, description="Minumum p-Phase current limit (pu)"]
        I_pre(t), [guess=0.8574182321, description=""]
        I_post(t), [guess=0.9259688073, description=""]
        VDL1_out(t), [guess=1.1, description=""]
        VDL2_out(t), [guess=1.2, description=""]
    end
    @equations begin
        Voltage_dip ~ ifelse(Vt_in.u<V_dip, 1, ifelse(Vt_in.u>V_up, 1, 0))

        #T_rv * Dt(V_tfilt) ~ Vt_in.u - V_tfilt
        simpleLag.in ~ Vt_in.u
        V_tfilt ~ simpleLag.out

        V_tfiltlim ~ max(V_tfilt, 0.01) #lowlimit(V_tfilt, 0.01)
        #q-phase current
        ΔV_t ~ V_ref0 - V_tfilt

        #ΔV_tdbd ~ deadband(ΔV_t, dbd1, dbd2)
        deadband.in ~ ΔV_t
        ΔV_tdbd ~ deadband.out

        I_qinj ~ (1-Voltage_dip) * clamp(K_qv*ΔV_tdbd, I_ql1, I_qh1) #limiter(K_qv*ΔV_tdbd, I_ql1, I_qh1)

        #T_p * Dt(P_PF) ~ P_e.u - P_PF
        simpleLag1.in ~ P_e.u
        P_PF ~ simpleLag1.out

        Q_con ~ ifelse(PfFlag, P_PF * tan(P_faref.u), Qext_in.u) #PfFlag * P_PF * tan(P_faref.u) + (1-PfFlag) * Qext_in.u
        Q_lim ~ clamp(Q_con, Q_min, Q_max) #limiter(Q_con, Q_min, Q_max)
        ΔQ ~ Q_lim - Q_gen.u
        Dt(s_Q) ~ K_qi * (1-Voltage_dip)* ΔQ
        V_in ~ K_qp * ΔQ + s_Q
        V_lima ~ clamp(V_in, V_min, V_max) #limiter(V_in, V_min, V_max)
        V_mod ~ ifelse((!PfFlag) && (!Vflag) && QFlag, (V_0 - PfFlag), V_bias) #ifelse((1-PfFlag) * (1-Vflag) * QFlag > 0, (V_0 - PfFlag), V_bias)
        V_con ~ ifelse(Vflag, V_lima, (Q_con + V_mod)) #Vflag * V_lima + (1-Vflag) * (Q_con + V_mod)
        V_limb ~ clamp(V_con, V_min, V_max) #limiter(V_con, V_min, V_max)
        ΔV ~ V_limb - V_tfilt
        Dt(s_V) ~ (1-Voltage_dip) * K_vi * ΔV
        I_in ~ K_vp * ΔV + s_V
        I_lim ~ clamp(I_in, I_qmin, I_qmax) #limiter(I_in, I_qmin, I_qmax)
        I_t ~ Q_con / V_tfiltlim
        ΔI ~ I_t - I_qin
        T_iq * Dt(I_qin) ~ ΔI
        I_qcon ~ ifelse(QFlag, I_lim, I_qin) #QFlag * I_lim + (1-QFlag) * I_qin
        I_sum ~ I_qcon + I_qinj
        I_qcmd ~ clamp(I_sum, I_qmin, I_qmax) #limiter(I_sum, I_qmin, I_qmax)
        #p-phase current
        P_in ~ ifelse(PfFlag, Wg.u * Pref_in.u, Pref_in.u) #(PfFlag * Wg.u + 1-PfFlag) * Pref_in.u
        ΔP ~ P_in - P_refout
        ΔP_lim ~ clamp(ΔP, dP_min, dP_max) #limiter(ΔP, dP_min, dP_max)
        T_pord * Dt(P_refout) ~ ΔP_lim
        P_lim ~ clamp(P_refout, P_min, P_max) #limiter(P_refout, P_min, P_max)
        I_pref ~ P_lim/V_tfiltlim
        I_pcmd ~ clamp(I_pref, I_pmin, I_pmax) #limiter(I_pref, I_pmin, I_pmax)
        #VDL tables
        VDL1_out ~ VDL(V_tfilt, Vq1, Vq2, Vq3, Vq4, Iq1, Iq2, Iq3, Iq4)
        VDL2_out ~ VDL(V_tfilt, Vp1, Vp2, Vp3, Vp4, Ip1, Ip2, Ip3, Ip4)
        #current limiter logic
        I_pmin ~ 0
        I_qmin ~ -I_qmax
        I_pre ~ ifelse(PqFlag, (sqrt(I_max) - sqrt(abs(I_pcmd))), (sqrt(I_max) - sqrt(abs(I_qcmd)))) #(1-PqFlag) * (sqrt(I_max) - sqrt(abs(I_qcmd))) + PqFlag * (sqrt(I_max) - sqrt(abs(I_pcmd)))
        I_post ~ ifelse(I_pre<0, 0, sqrt(I_pre))
        I_pmax ~ ifelse(PqFlag, min(VDL2_out, I_max), min(VDL2_out, I_post)) #PqFlag * min(VDL2_out, I_max) + (1-PqFlag) * min(VDL2_out, I_post)
        I_qmax ~ ifelse(PqFlag, min(VDL1_out, I_post), min(VDL1_out, I_max)) #PqFlag * min(VDL1_out, I_post) + (1-PqFlag) * min(VDL1_out, I_max)
        #outputs
        Iqcmd_out.u ~ I_qcmd
        Ipcmd_out.u ~ I_pcmd
    end
end


@mtkmodel reec_b_pf begin
    @structural_parameters begin
        PfFlag = false
        Vflag = false
        QFlag = false
        PqFlag = false
    end
    @parameters begin
        V_dip, [description="Low voltage condition trigger voltage (pu)"]
        V_up, [description="High voltage condition trigger voltage (pu)"]
        T_rv, [description="Terminal bus voltage filter time constant (s)"]
        V_ref0, [description="Reference voltage for reactive current injection (pu)"]
        dbd1, [description="Overvoltage deadband for reactive current injection (pu)"]
        dbd2, [description="Undervoltage deadband for reactive current injection (pu)"]
        K_qv, [description="Reactive current injection gain (pu/pu)"]
        I_qh1, [description="Maximum reactive current injection (pu on mbase)"]
        I_ql1, [description="Minimum reactive current injection (pu on mbase)"]
        T_p, [description="Active power filter time constant (s)"]
        Q_min, [description="Minimum reactive power when Vflag = 1 (pu on mbase)"]
        Q_max, [description="Maximum reactive power when Vflag = 1 (pu on mbase)"]
        V_min, [description="Minimum voltage at inverter terminal bus (pu)"]
        V_max, [description="Maximum voltage at inverter terminal bus (pu)"]
        K_qp, [description="Local Q regulator proportional gain (pu/pu)"]
        K_qi, [description="Local Q regulator integral gain (pu/pu-s)"]
        K_vp, [description="Local voltage regulator proportional gain (pu/pu)"]
        K_vi, [description="Local voltage regulator integral gain (pu/pu-s)"]
        I_max, [description="Maximum apparent current (pu on mbase)"]
        T_iq, [description="Reactive current regulator lag time constant (s)"]
        T_pord, [description="Inverter power order lag time constant (s)"]
        P_min, [description="Minimum active power (pu on mbase)"]
        P_max, [description="Maximum active power (pu on mbase)"]
        dP_min, [description="Active power down-ramp limit (pu/s on mbase)"]
        dP_max, [description="Active power up-ramp limit (pu/s on mbase)"]
    end
    @components begin
        Vt_in = RealInput(guess=1.001047)
        P_e = RealInput(guess=0.800721) #Inverter active power (pu on mbase)
        P_initial = RealInput(guess=0.800721)
        Q_initial = RealInput(guess=-0.30027)
        Qext_in = RealInput(guess=-0.30027)
        Pref_in = RealInput(guess=0.800721)
        Q_gen = RealInput(guess=-0.30027)
        # outputs
        Iqcmd_out = RealOutput(guess=-0.299956)
        Ipcmd_out = RealOutput(guess=0.799884)
        Vt_filt = RealOutput(guess=1.001047)

        #building blocks
        simpleLag = PowerDynamics.Library.SimpleLag(K=1, T=T_rv, guess=1.001047)
        simpleLag1 = PowerDynamics.Library.SimpleLag(K=1, T=T_p, guess=0.800721)
        deadband = PowerDynamics.Library.DeadZone(uMax=dbd2, uMin=dbd1)
        PI_freeze = PowerDynamics.Library.P_I_Lim_freeze(K_p=K_qp, K_i=K_qi, T=1, outMin=V_min, outMax=V_max, guess=1.001047, guessx=1.001047)
        PI_freeze_var = PowerDynamics.Library.P_I_varLim_freeze(K_p=K_vp, K_i=K_vi, T=1, guess=-0.299956, guessx=-0.299956)
        simpleLag_freeze = PowerDynamics.Library.SimpleLag_freeze(K=1, T=T_iq, guess=-0.299956, guessx=-0.299956)
        P_limLag = PowerDynamics.Library.SimpleLag_2Lims_freeze(K=1, T=T_pord, doutMin=dP_min, doutMax=dP_max, outMin=P_min, outMax=P_max, guess=0.800721, guessin=0.800721, guessx=0.800721)
    end
    @variables begin
        Voltage_dip(t), [guess=0, description="freeze states if Voltagedip=1"]
        V_tfilt(t), [guess=1.001047, description="Voltage after filter"]
        V_tfiltlim(t), [guess=1.001047, description="Voltage after filter with lower limit 0.01"]
        ΔV_t(t), [guess=1.001047, description="Difference between filterd terminal voltage and reference voltage"]
        ΔV_tdbd(t), [guess=-0.981046, description="Voltage after deadband"]
        I_qinj(t), [guess=0, description="Limited current injection q-Phase from Voltage"]
        P_PF(t), [guess=0.800721, description="Inverter active power after filter"]
        Q_con(t), [guess=-0.30027, description="Reactive Power after PfFlag"]
        Q_lim(t), [guess=-0.30027, description="ReactivePower after limiter"]
        ΔQ(t), [guess=0, description="Difference between Q_lim and Q_gen"]
        #s_Q(t), [guess=4.8269116e-9, description="Frozen state in Q regulator"]
        #s_Qint(t), [guess=1, description=""]
        #V_in(t), [guess=1, description="Voltage after local Q regulator"]
        V_lima(t), [guess=1.001047, description="Limited voltage after Q regulator"]
        V_con(t), [guess=1.001047, description="Voltage after Vflag"]
        V_limb(t), [guess=1.001047, description="Limited voltage V_con"]
        ΔV(t), [guess=0, description="Difference between V_limb and V_tfilt"]
        #s_V(t), [guess=6.6613381e-14, description="Frozen state at local voltage regulator"]
        #s_Vint(t), [guess=-0.0567, description=""]
        #I_in(t), [guess=-0.0567, description="Current after local voltage regulator"]
        I_lim(t), [guess=-0.299956, description="limited current after voltage regulator"]
        I_t(t), [guess=-0.299956, description="Current from Q_con/V_tfiltlim"]
        #ΔI(t), [guess=-1.9310942e-14, description=""]
        I_qin(t), [guess=-0.299956, description="Current after Reactive current regulator"]
        I_qcon(t), [guess=-0.299956, description="Current after QFlag"]
        I_sum(t), [guess=-0.299956, description="sum of I_qcon and I_qinj"]
        I_qcmd(t), [guess=-0.299956, description="q-Phase output current"]
        #P_refout(t), [guess=0.015, description="Active power after inverter power order"]
        P_lim(t), [guess=0.800721, description="Limited active power after inverter power order"]
        #ΔP(t), [guess=-7.3598558e-12, description="Active power difference between P_ref and P_refout"]
        #ΔP_lim(t), [guess=-7.3598558e-12, description="Ramp-limited active power difference"]
        I_pref(t), [guess=0.799884, description="Current from P_lim/V_tfiltlim"]
        I_pcmd(t), [guess=0.799884, description="p-Phase output current"]
        I_qmin(t), [guess=-1.25, description="Minumum q-Phase current limit (pu)"]
        I_qmax(t), [guess=1.25, description="Maximum q-Phase current limit (pu)"]
        I_pmax(t), [guess=1.213477, description="Maximum p-Phase current limit (pu)"]
        I_pmin(t), [guess=0, description="Minumum p-Phase current limit (pu)"]
    end
    @equations begin
        #T_rv * Dt(V_tfilt) ~ Vt_in.u - V_tfilt
        simpleLag.in ~ Vt_in.u
        V_tfilt ~ simpleLag.out
        Voltage_dip ~ ifelse(V_tfilt<V_dip, 1, ifelse(V_tfilt>V_up, 1, 0)) #TODO in PF wird Voltage_dip mit V_tfilt bestimmt

        V_tfiltlim ~ max(V_tfilt, 0.01) #lowlimit(V_tfilt, 0.01)
        #q-phase current
        ΔV_t ~ - V_ref0 + V_tfilt #TODO in PF VZ anders im Vgl zur angegebenen Quelle

        #ΔV_tdbd ~ deadband(ΔV_t, dbd1, dbd2)
        deadband.in ~ ΔV_t
        ΔV_tdbd ~ deadband.out

        I_qinj ~ clamp(K_qv*ΔV_tdbd, I_ql1, I_qh1) #limiter(K_qv*ΔV_tdbd, I_ql1, I_qh1)

        #T_p * Dt(P_PF) ~ P_e.u - P_PF
        simpleLag1.in ~ P_e.u
        P_PF ~ simpleLag1.out

        Q_con ~ ifelse(PfFlag, P_PF * Q_initial.u/P_initial.u, Qext_in.u)  #ifelse(PfFlag, P_PF * tan(P_faref.u), Qext_in.u)
        Q_lim ~ clamp(Q_con, Q_min, Q_max) #limiter(Q_con, Q_min, Q_max)
        ΔQ ~ Q_lim - Q_gen.u
        PI_freeze.in ~ ΔQ
        PI_freeze.freeze ~ Voltage_dip
        V_lima ~ PI_freeze.out

        V_con ~ ifelse(Vflag, V_lima, Q_con) #Vflag * V_lima + (1-Vflag) * V_ref0
        V_limb ~ clamp(V_con, V_min, V_max) #limiter(V_con, V_min, V_max)
        ΔV ~ V_limb - V_tfilt
        PI_freeze_var.in ~ ΔV
        PI_freeze_var.freeze ~ Voltage_dip
        PI_freeze_var.outMin ~ I_qmin
        PI_freeze_var.outMax ~ I_qmax
        I_lim ~ PI_freeze_var.out

        I_t ~ Q_con / V_tfiltlim
        simpleLag_freeze.in ~ I_t
        simpleLag_freeze.freeze ~ Voltage_dip
        I_qin ~ simpleLag_freeze.out
        I_qcon ~ ifelse(QFlag, I_lim, I_qin)

        I_sum ~ I_qcon + I_qinj
        I_qcmd ~ limiter(I_sum, I_qmin, I_qmax) #clamp(I_sum, I_qmin, I_qmax) #limiter(I_sum, I_qmin, I_qmax)

        #p-phase current
        P_limLag.in ~ Pref_in.u
        P_limLag.freeze ~ Voltage_dip
        P_lim ~ P_limLag.out
        I_pref ~ P_lim/V_tfiltlim
        I_pcmd ~ limiter(I_pref, I_pmin, I_pmax) #clamp(I_pref, I_pmin, I_pmax) #limiter(I_pref, I_pmin, I_pmax)

        #current limiter logic
        I_pmin ~ 0
        I_qmin ~ - I_qmax
        I_pmax ~ ifelse(PqFlag, I_max, sqrt(I_max^2 - min(I_qcmd^2, I_max^2)))
        I_qmax ~ ifelse(PqFlag, sqrt(I_max^2 - min(I_pcmd^2, I_max^2)), I_max)
        #outputs
        Iqcmd_out.u ~ I_qcmd
        Ipcmd_out.u ~ I_pcmd
        Vt_filt.u ~ V_tfilt
    end
end
