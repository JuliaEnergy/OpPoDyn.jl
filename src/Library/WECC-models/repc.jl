@mtkmodel repc_a begin
    @structural_parameters begin
        RefFlag=false
        VcombFlag=false
        freqFlag=false
    end
    @components begin
        V_reg = RealInput(guess=1) #from bus - was ist Unterschied zu Vt = raw terminal voltage ? initialized based on initial power flow solution for V_reg and Q_gen
        Q_branch = RealInput(guess=-0.0567) #reactive power through a defined branch /from aggregated turbine model or collection point of wind plant
        P_branch = RealInput(guess=0.015) #from aggregated turbine model or collection point of wind plant
        freq = RealInput(guess=60) #from aggregated turbine model or collection point of wind plant
        V_ref = RealInput(guess=1)
        Q_ref = RealInput(guess=-0.0567)
        V_diff = RealInput(guess=1)
        # outputs
        Pref_out = RealOutput() # active power [pu]
        Qext_out = RealOutput() # reactive power [pu]
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
    @variables begin
        Voltage_dip(t), [guess=0, description="freeze states if Voltagedip=1"]
        V_droop(t), [guess=0.998, description="Input voltage for VcombFlag=0"]
        V_in(t), [guess=1, description="resulting voltage after VcombFlag"]
        V_fltr(t), [guess=0.99887, description="Voltage after filter"]
        ΔV(t), [guess=0, description="Voltage difference between filter and reference"]
        Q_fltr(t), [guess=-0.056657, description="Reactive power after filter"]
        ΔQ(t), [guess=0, description="Reactive power differemce between filter and reference"]
        ΔQ_in(t), [guess=0, description="Input depending on RefFlag=1 (ΔV) or RefFlag=0 (ΔQ)"]
        ΔQ_dbd(t), [guess=0, description="Reactive power differece after deadband"]
        Q_e(t), [guess=0, description="Reactive power after error limits"]
        Q_x(t), [guess=-0.0567, description="Reactive power before limits"]
        Q_res(t), [guess=0, description=""]
        Q_I(t), [guess=-0.056635, description=""]
        Q_lim(t), [guess=-0.0567, description="Reactive power after reactive power limits"]
        Q_ext(t), [guess=-0.056657, description="Reactive power output"]
        Δf_deadband(t), [guess=0, description="frequency difference after deadband"]
        Δf_corr(t), [guess=0, description="Frequency difference after droop"]
        P_branchp(t), [guess=0.015, description="Active power after T_p"]
        f_e(t), [guess=0, description="frequency after frequency limits"]
        P_e(t), [guess=0.015,description="Active power before power limits"]
        P_lim(t), [guess=0015, description="Active power after power limits"]
        P_refa(t), [guess=0.015, description="Active Power reference if Freq_flag=1"]
        P_ref(t), [guess=0.015, description="Active power output"]
    end
    @equations begin
        Voltage_dip ~ ifelse(V_reg.u<V_frz, 1,0) #0
        #reactive power control
        V_droop ~ K_c * Q_branch.u + V_reg.u
        V_in ~ ifelse(VcombFlag, V_diff.u, V_droop) #VcombFlag * V_diff.u + (1-VcombFlag) * V_droop
        T_fltr * Dt(V_fltr) ~ V_in - V_fltr
        ΔV ~ V_ref.u - V_fltr
        T_fltr * Dt(Q_fltr) ~ Q_branch.u - Q_fltr
        ΔQ ~ Q_ref.u - Q_fltr
        ΔQ_in ~ ifelse(RefFlag, ΔV, ΔQ) #RefFlag * ΔV + (1-RefFlag)*ΔQ
        ΔQ_dbd ~ deadband(ΔQ_in, dbd_dn, dbd_up)
        Q_e ~ limiter(ΔQ_dbd, e_min, e_max)
        Q_res ~ (1-Voltage_dip) * Q_e
        Dt(Q_I) ~ K_i * Q_res
        Q_x ~ K_p * Q_e + Q_I
        Q_lim ~ limiter(Q_x, Q_min, Q_max)
        T_fv * Dt(Q_ext) - T_ft* Dt(Q_lim) ~ Q_lim - Q_ext
        #active power control
        Δf_deadband ~ deadband(freq_ref-freq.u, fdbd1, fdbd2)
        Δf_corr ~ lowlimit(Δf_deadband * D_up, 0) + uplimit(Δf_deadband * D_dn, 0)
        T_p * Dt(P_branchp) ~ P_branch.u - P_branchp
        f_e ~ limiter(P_plantref-P_branchp+Δf_corr, femin, femax)
        Dt(P_e) - Dt(f_e) * K_pg ~ K_ig * f_e
        P_lim ~ limiter(P_e, P_min, P_max)
        T_lag * Dt(P_refa) ~ P_lim - P_refa
        P_ref ~ ifelse(freqFlag, P_refa, p_0) #freqFlag * P_refa + (1-freqFlag) * p_0
        #outputs
        Pref_out.u ~ P_ref
        Qext_out.u ~ Q_ext
    end
end


