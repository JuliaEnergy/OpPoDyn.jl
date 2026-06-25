@mtkmodel repc_a begin
    @structural_parameters begin
        RefFlag=false
        VcombFlag=false
        freqFlag=false
    end
    @parameters begin
        K_p, [description="Proportional gain in pu/pu"]
        K_i, [description="Integral gain in pu/pu "]
        T_fltr, [description="Voltage or reactive power measurement filter time constant in s"]
        T_ft, [description="Lead time constant in s"]
        T_fv, [description="Lag time constant in s"]
        V_frz, [description="Voltage below which plant control integrator state (s2) is frozen in pu"]
        R_c, [description="Line drop compensation resistance in pu"]
        X_c, [description="Current compensation constant (to emulate droop or line drop compensation) in pu"]
        K_c, [description="Gain on reactive current compensation in pu"]
        e_max, [description="Maximum error limit in pu"]
        e_min, [description="Minimum error limit in pu"]
        dbd_up, [description="Deadband in control in pu"]
        dbd_dn, [description="Deadband in control in pu"]
        Q_max, [description="Maximum Q control output in pu"]
        Q_min, [description="Minimum Q control output in pu"]
        K_pg, [description="Proportional gain for power control in pu/pu"]
        K_ig, [description="Integral gain for power control in pu/pu"]
        T_p, [description="Lag time constant on Pgen measurement in s"]
        fdbd1, [description="Deadband downside in pu"]
        fdbd2, [description="Deadband upside in pu"]
        femax, [description="Maximum error limit in pu"]
        femin, [description="Minimum error limit in pu"]
        P_max, [description="Maximum Power in pu"]
        P_min, [description="Minimum Power in pu"]
        T_lag, [description="Lag time constant in s"]
        D_dn, [description="Downside droop in pu/pu"]
        D_up, [description="upside droop in pu/pu"]
        P_plantref, [description="Initial power reference in pu"]
        freq_ref, [description="Frequency reference"]
        p_0, [guess=0.015, description="Initial active power in pu"]
    end
    @components begin
        V_reg = RealInput(guess=1)
        Q_branch = RealInput(guess=-0.056656801)
        P_branch = RealInput(guess=0.015)
        freq = RealInput(guess=60)
        V_ref = RealInput(guess=1)
        Q_ref = RealInput(guess=-0.056658)
        V_diff = RealInput(guess=1.0001042)
        # outputs
        Pref_out = RealOutput(guess=0.015)
        Qext_out = RealOutput(guess=-0.056656797)

        #building blocks
        simpleLag2 = PowerDynamics.Library.SimpleLag(K=1, T=T_fltr, guess=-0.056656813)
        simpleLag3 = PowerDynamics.Library.SimpleLag(K=1, T=T_lag, guess=0.99886686)
        leadLag = PowerDynamics.Library.LeadLag(K=1, T1=T_ft, T2=T_fv, guess=-0.056656797)
        deadband = PowerDynamics.Library.DeadZone(uMax=dbd_up, uMin=dbd_dn)
        if freqFlag
            simpleLag = PowerDynamics.Library.SimpleLag(K=1, T=T_p, guess=0.015)
            simpleLag1 = PowerDynamics.Library.SimpleLag(K=1, T=T_lag, guess=0.015)
            f_deadband = PowerDynamics.Library.DeadZone(uMax=fdbd2, uMin=fdbd1)
        end
    end
    @variables begin
        Voltage_dip(t), [guess=0, description="freeze states if Voltagedip=1"]
        V_droop(t), [guess=0.99886686, description="Input voltage for VcombFlag=0"]
        V_in(t), [guess=0.99886686, description="resulting voltage after VcombFlag"]
        V_fltr(t), [guess=0.99886686, description="Voltage after filter"]
        ΔV(t), [guess=0.001133136, description="Voltage difference between filter and reference"]
        Q_fltr(t), [guess=-0.056656813, description="Reactive power after filter"]
        ΔQ(t), [guess=-1.1867191e-6, description="Reactive power differemce between filter and reference"]
        ΔQ_in(t), [guess=-1.1867191e-6, description="Input depending on RefFlag=1 (ΔV) or RefFlag=0 (ΔQ)"]
        ΔQ_dbd(t), [guess=-1.1867191e-6, description="Reactive power differece after deadband"]
        Q_e(t), [guess=-1.1867191e-6, description="Reactive power after error limits"]
        Q_x(t), [guess=-0.056656797, description="Reactive power before limits"]
        Q_res(t), [guess=-1.1867191e-6, description=""]
        Q_I(t), [guess=-0.056635436, description=""]
        Q_lim(t), [guess=-0.056656797, description="Reactive power after reactive power limits"]
        Q_ext(t), [guess=-0.056656797, description="Reactive power output"]
        if freqFlag
            Δf_deadband(t), [guess=0, description="frequency difference after deadband"]
            Δf_corr(t), [guess=0, description="Frequency difference after droop"]
            P_branchp(t), [guess=0.015, description="Active power after T_p"]
            f_e(t), [guess=1.4455451e-14, description="frequency after frequency limits"]
            P_e(t), [guess=0.015, description="Active power before power limits"]
            P_lim(t), [guess=0.015, description="Active power after power limits"]
            P_refa(t), [guess=0.015, description="Active Power reference if Freq_flag=1"]
        end
        P_ref(t), [guess=0.015, description="Active power output"]
    end
    @equations begin
        Voltage_dip ~ ifelse(V_reg.u<V_frz, 1,0)
        #reactive power control
        V_droop ~ K_c * Q_branch.u + V_reg.u
        if VcombFlag
            V_in ~ V_diff.u
        else
            V_in ~ V_droop
        end

        simpleLag3.in ~ V_in
        V_fltr ~ simpleLag3.out

        ΔV ~ V_ref.u - V_fltr

        simpleLag2.in ~ Q_branch.u
        Q_fltr ~ simpleLag2.out

        ΔQ ~ Q_ref.u - Q_fltr
        if RefFlag
            ΔQ_in ~ ΔV
        else
            ΔQ_in ~ ΔQ
        end

        deadband.in ~ ΔQ_in
        ΔQ_dbd ~ deadband.out

        Q_e ~ clamp(ΔQ_dbd, e_min, e_max)
        Q_res ~ (1-Voltage_dip) * Q_e
        Dt(Q_I) ~ K_i * Q_res
        Q_x ~ K_p * Q_e + Q_I
        Q_lim ~ clamp(Q_x, Q_min, Q_max)

        leadLag.in ~ Q_lim
        Q_ext ~ leadLag.out

        #active power control
        if freqFlag
            f_deadband.in ~ freq_ref-freq.u
            Δf_deadband ~ f_deadband.out
            Δf_corr ~ max(Δf_deadband * D_up, 0) + min(Δf_deadband * D_dn, 0)
            simpleLag.in ~ P_branch.u
            P_branchp ~ simpleLag.out
            f_e ~ clamp(P_plantref-P_branchp+Δf_corr, femin, femax)
            Dt(P_e) - Dt(f_e) * K_pg ~ K_ig * f_e
            P_lim ~ clamp(P_e, P_min, P_max)
            simpleLag1.in ~ P_lim
            P_refa ~ simpleLag1.out
            P_ref ~ P_refa
        else
            P_ref ~ p_0
        end
        #outputs
        Pref_out.u ~ P_ref
        Qext_out.u ~ Q_ext
    end
end


@mtkmodel repc_a_pf begin
    @structural_parameters begin
        RefFlag=false
        VcombFlag=false
        freqFlag=false
    end
    @parameters begin
        K_p, [description="Proportional gain in pu/pu"]
        K_i, [description="Integral gain in pu/pu "]
        T_fltr1, [description="Voltage or reactive power measurement filter time constant in s"]
        T_fltr2, [description="Voltage or reactive power measurement filter time constant in s"]
        T_ft, [description="Lead time constant in s"]
        T_fv, [description="Lag time constant in s"]
        V_frz, [description="Voltage below which plant control integrator state (s2) is frozen in pu"]
        R_c, [description="Line drop compensation resistance in pu"]
        X_c, [description="Current compensation constant (to emulate droop or line drop compensation) in pu"]
        K_c, [description="Gain on reactive current compensation in pu"]
        e_max, [description="Maximum error limit in pu"]
        e_min, [description="Minimum error limit in pu"]
        dbd_up, [description="Deadband in control in pu"]
        dbd_dn, [description="Deadband in control in pu"]
        Q_max, [description="Maximum Q control output in pu"]
        Q_min, [description="Minimum Q control output in pu"]
        K_pg, [description="Proportional gain for power control in pu/pu"]
        K_ig, [description="Integral gain for power control in pu/pu"]
        T_p, [description="Lag time constant on Pgen measurement in s"]
        fdbd1, [description="Deadband downside in pu"]
        fdbd2, [description="Deadband upside in pu"]
        femax, [description="Maximum error limit in pu"]
        femin, [description="Minimum error limit in pu"]
        P_max, [description="Maximum Power in pu"]
        P_min, [description="Minimum Power in pu"]
        T_lag, [description="Lag time constant in s"]
        D_dn, [description="Downside droop in pu/pu"]
        D_up, [description="upside droop in pu/pu"]
        #P_plantref, [description="Initial power reference in pu"]
        freq_ref, [description="Frequency reference"]
    end
    @components begin
        V_reg_re = RealInput(guess=1.001044) #from bus - was ist Unterschied zu Vt = raw terminal voltage ? initialized based on initial power flow solution for V_reg and Q_gen
        V_reg_im = RealInput(guess=0.00231)
        P_plantref = RealInput(guess=0.800721) #TODO Wert in PF
        Q_branch = RealInput(guess=-0.30027) #reactive power through a defined branch /from aggregated turbine model or collection point of wind plant
        P_branch = RealInput(guess=0.800721) #from aggregated turbine model or collection point of wind plant
        I_branch = RealInput(guess=0.854276)
        freq = RealInput(guess=1) #from aggregated turbine model or collection point of wind plant
        V_ref = RealInput(guess=0.903563)
        Q_ref = RealInput(guess=-0.30027)
        # outputs
        Pref_out = RealOutput(guess=0.800721) # active power [pu]
        Qext_out = RealOutput(guess=-0.30027) # reactive power [pu]

        #building blocks
        simpleLag2 = PowerDynamics.Library.SimpleLag(K=1, T=T_fltr2, guess=-0.30027)
        simpleLag3 = PowerDynamics.Library.SimpleLag(K=1, T=T_fltr1, guess=1.001047)
        leadLag = PowerDynamics.Library.LeadLag(K=1, T1=T_ft, T2=T_fv, guess=-0.30027)
        deadband = PowerDynamics.Library.DeadZone(uMax=dbd_up, uMin=dbd_dn)
        PI_lim_Q = PowerDynamics.Library.P_I_Lim_freeze(K_p=K_p, K_i=K_i, T=1, outMin=Q_min, outMax=Q_max, guess=-0.30027, guessx=-0.30027)
        if freqFlag
            simpleLag = PowerDynamics.Library.SimpleLag(K=1, T=T_p, guess=0.800721)
            simpleLag1 = PowerDynamics.Library.SimpleLag(K=1, T=T_lag, guess=0.800721)
            f_deadband = PowerDynamics.Library.DeadZone(uMax=fdbd2, uMin=fdbd1)
            PI_lim_P = PowerDynamics.Library.P_I_Lim(K_p=K_pg, K_i=K_ig, T=1, outMin=P_min, outMax=P_max, guess=0.800721, guessx=0.800721)
        end
    end
    @variables begin
        Voltage_dip(t), [guess=0, description="freeze states if Voltagedip=1"]
        V_reg(t), [guess=1.001047, description="absolute value of V_reg"]
        V_diff(t), [guess=1.001047, description=""]
        V_droop(t), [guess=0.995041, description="Input voltage for VcombFlag=0"]
        V_in(t), [guess=1.001047, description="resulting voltage after VcombFlag"]
        V_fltr(t), [guess=1.001047, description="Voltage after filter"]
        ΔV(t), [guess=0, description="Voltage difference between filter and reference"]
        Q_fltr(t), [guess=-0.30027, description="Reactive power after filter"]
        ΔQ(t), [guess=0, description="Reactive power differemce between filter and reference"]
        ΔQ_in(t), [guess=0, description="Input depending on RefFlag=1 (ΔV) or RefFlag=0 (ΔQ)"]
        ΔQ_dbd(t), [guess=0, description="Reactive power differece after deadband"]
        Q_e(t), [guess=0, description="Reactive power after error limits"]
        #Q_x(t), [guess=-0.056656797, description="Reactive power before limits"]
        #Q_res(t), [guess=-1.1867191e-6, description=""]
        #Q_I(t), [guess=-0.056635436, description=""]
        Q_lim(t), [guess=-0.30027, description="Reactive power after reactive power limits"]
        Q_ext(t), [guess=-0.30027, description="Reactive power output"]
        if freqFlag
            Δf_deadband(t), [guess=0, description="frequency difference after deadband"]
            Δf_corr(t), [guess=0, description="Frequency difference after droop"]
            P_branchp(t), [guess=0.800721, description="Active power after T_p"]
            f_e(t), [guess=0, description="frequency after frequency limits"]
            P_lim(t), [guess=0.800721, description="Active power after power limits"]
            P_refa(t), [guess=0.800721, description="Active Power reference if Freq_flag=1"]
        end
        P_ref(t), [guess=0.800721, description="Active power output"]
    end
    @equations begin
        V_reg ~ sqrt(V_reg_re.u^2 + V_reg_im.u^2)
        Voltage_dip ~ ifelse(V_reg<V_frz, 1,0) #0

        #reactive power control
        V_diff ~ sqrt((V_reg_re.u-(R_c*I_branch.u))^2 + (V_reg_im.u-(X_c*I_branch.u))^2)
        V_droop ~ K_c * Q_branch.u + V_reg
        V_in ~ ifelse(VcombFlag, V_diff, V_droop)

        #T_fltr * Dt(V_fltr) ~ V_in - V_fltr
        simpleLag3.in ~ V_in
        V_fltr ~ simpleLag3.out

        ΔV ~ V_ref.u - V_fltr

        #T_fltr * Dt(Q_fltr) ~ Q_branch.u - Q_fltr
        simpleLag2.in ~ Q_branch.u
        Q_fltr ~ simpleLag2.out

        ΔQ ~ Q_ref.u - Q_fltr
        ΔQ_in ~ ifelse(RefFlag, ΔV, ΔQ) #RefFlag * ΔV + (1-RefFlag)*ΔQ

        #ΔQ_dbd ~ deadband(ΔQ_in, dbd_dn, dbd_up)
        deadband.in ~ ΔQ_in
        ΔQ_dbd ~ deadband.out

        Q_e ~ clamp(ΔQ_dbd, e_min, e_max) #limiter(ΔQ_dbd, e_min, e_max)
        PI_lim_Q.in ~ Q_e
        PI_lim_Q.freeze ~ Voltage_dip
        Q_lim ~ PI_lim_Q.out

        #T_fv * Dt(Q_ext) - T_ft* Dt(Q_lim) ~ Q_lim - Q_ext
        leadLag.in ~ Q_lim
        Q_ext ~ leadLag.out

        #active power control
        if freqFlag
            f_deadband.in ~ freq_ref-freq.u
            Δf_deadband ~ f_deadband.out
            Δf_corr ~ max(Δf_deadband * D_up, 0) + min(Δf_deadband * D_dn, 0)
            simpleLag.in ~ P_branch.u
            P_branchp ~ simpleLag.out
            f_e ~ clamp(P_plantref.u-P_branchp+Δf_corr, femin, femax)
            PI_lim_P.in ~ f_e
            P_lim ~ PI_lim_P.out
            simpleLag1.in ~ P_lim
            P_refa ~ simpleLag1.out
            P_ref ~ P_refa
        else
            P_ref ~ P_plantref.u
        end
        #outputs
        Pref_out.u ~ P_ref
        Qext_out.u ~ Q_ext
    end
end


