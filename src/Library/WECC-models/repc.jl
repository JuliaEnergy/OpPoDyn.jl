@mtkmodel repc_a begin
    @structural_parameters begin
        #vf_input = true
    end
    @components begin
        #terminal=Terminal()
        # inputs
        #if vf_input
        #    vf_in = RealInput(guess=1) # field voltage input [pu]
        #end
        V_reg = RealInput(guess=1) #from bus - was ist Unterschied zu Vt = raw terminal voltage ? initialized based on initial power flow solution for V_reg and Q_gen
        I_branch = RealInput(guess=1) #current through a defined branch
        Q_branch = RealInput(guess=1) #reactive power through a defined branch /from aggregated turbine model or collection point of wind plant
        P_branch = RealInput(guess=1) #from aggregated turbine model or collection point of wind plant
        freq = RealInput(guess=1) #from aggregated turbine model or collection point of wind plant
        V_ref = RealInput(guess=1)
        Q_ref = RealInput(guess=1)
        # outputs
        Pref_out = RealOutput() # active power [pu]
        Qext_out = RealOutput() # reactive power [pu]
    end
    @parameters begin
        #primary parameters
        K_p, [description="Proportional gain in pu/pu"]
        K_i, [description="Integral gain in pu/pu "]
        T_fltr, [description="Voltage or reactive power measurement filter time constant in s"]
        T_ft, [description="Lead time constant in s"]
        T_fv, [description="Lag time constant in s"]
        RefFlag=1, [description="1 for voltage control or 0 for reactive power control"]
        V_frz, [description="Voltage below which plant control integrator state (s2) is frozen in pu"]
        R_c, [description="Line drop compensation resistance in pu"]
        X_c, [description="Current compensation constant (to emulate droop or line drop compensation) in pu"]
        K_c, [description="Gain on reactive current compensation in pu"]
        VcombFlag, [description="Selection of droop (0) or line drop compensation (1)"]
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
        freq_ref, [description="Frequency reference in pu"]
        freqFlag, [description="Flag to turn on (1) or off (0) the active power control loop within the plant controller"]
        #vbus, [description="The bus number in powerflow from which V_reg, freq is picked up (i.e. the voltage being regulated and frequency being controlled; it can be the terminal of the aggregated WTG model or the point of interconnection)"]
        #branch, [description="The branch (actual definition depends on software program) from which I_branch, Q_branch and P_branch is being measured"]
        # input/parameter switches
        #if !vf_input
        #    vf_set, [guess=1, bounds=(0,Inf), description="field voltage"]
        #end
    end
    @variables begin
        V_comp(t), [description="Input voltage for VcombFlag=1"]
        V_droop(t), [description="Input voltage for VcombFlag=0"]
        V_in(t), [description="resulting voltage after VcombFlag"]
        V_fltr(t), [description="Voltage after filter"]
        ΔV(t), [description="Voltage difference between filter and reference"]
        Q_fltr(t), [description="Reactive power after filter"]
        ΔQ(t), [description="Reactive power differemce between filter and reference"]
        ΔQ_in(t), [description="Input depending on RefFlag=1 (ΔV) or RefFlag=0 (ΔQ)"]
        ΔQ_dbd(t), [description="Reactive power differece after deadband"]
        Q_e(t), [description="Reactive power after error limits"]
        Q_x(t), [description="Reactive power before limits"]
        s_2(t), [description="state of reactive power that can be frozen depending on V_reg"]
        Q_lim(t), [description="Reactive power after reactive power limits"]
        Q_ext(t), [description="Reactive power output"]
        Δf_deadband(t), [description="frequency difference after deadband"]
        Δf_corr(t), [description="Frequency difference after droop"]
        P_branchp(t), [description="Active power after T_p"]
        Δf_in(t), [description="Difference between P_plantref and Δf_corr and P_branchp"]
        f_e(t), [description="frequency after frequency limits"]
        P_e(t), [description="Active power before power limits"]
        P_lim(t), [description="Active power after power limits"]
        P_refa(t), [description="Active Power reference if Freq_flag=1"]
        P_ref(t), [description="Active power output"]
    end
    @equations begin
        #reactive power control
        V_droop ~ K_c * Q_branch.u + V_reg.u
        V_comp ~ abs(V_reg.u-(R_c+im*X_c)*I_branch.u)
        V_in ~ VcombFlag * V_comp + (1-VcombFlag)*V_droop
        T_fltr * Dt(V_fltr) ~ V_in - V_fltr
        ΔV ~ V_ref.u - V_fltr
        T_fltr * Dt(Q_fltr) ~ Q_branch.u - Q_fltr
        ΔQ ~ Q_ref.u - Q_fltr
        ΔQ_in ~ RefFlag * ΔV + (1-RefFlag)*ΔQ
        ΔQ_dbd ~ deadband(ΔQ_in, dbd_dn, dbd_up)
        Q_e ~ limiter(ΔQ_dbd, e_min, e_max)
        Q_x ~ Q_e * K_p + s_2
        Dt(s_2) ~ ifelse(V_reg.u<V_frz, 0, K_i*Q_e)
        Q_lim ~ limiter(Q_x, Q_min, Q_max)
        T_fv * Dt(Q_ext) ~ T_ft* Dt(Q_lim) + Q_lim - Q_ext
        #active power control
        Δf_deadband ~ deadband(freq_ref-freq.u, fdbd1, fdbd2)
        Δf_corr ~ lowlimit(Δf_deadband * D_up, 0) + uplimit(Δf_deadband * D_dn, 0) #min(0, Δf_deadband*D_dn) + max(0, Δf_deadband*D_up)
        T_p * Dt(P_branchp) ~ P_branch.u - P_branchp
        f_e ~ limiter(P_plantref-P_branchp+Δf_corr, femin, femax)
        Dt(P_e) ~ Dt(f_e) * K_pg + K_ig * f_e
        P_lim ~ limiter(P_e, P_min, P_max)
        T_lag * Dt(P_refa) ~ P_lim - P_refa
        P_ref ~ freqFlag * P_refa + (1-freqFlag) * P_plantref #oder P_branch?
        #outputs
        Pref_out.u ~ P_ref
        Qext_out.u ~ Q_ext
        #v_mag_out.u ~ v_mag
        #δout.u ~ δ
        #ωout.u ~ ω
    end
end


