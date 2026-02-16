# WECC Model as used in https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7741608&tag=1

@mtkmodel WECC_large_PV_pf begin
    @components begin
        terminal=Terminal()
        regca = regc_a_pf(
            I_qrmax = 9999,
            I_qrmin = -9999,
            T_g     = 0.02,
            T_fltr  = 0.02,
            Brkpt   = 0.9,
            Zerox   = 0.5,
            L_vpl1  = 1.22,
            rrpwr   = 10,
            V_0lim = 1.2,
            K_hv = 0.7,
            lvpnt0 = 0.4,
            lvpnt1 = 0.8,
            I_0lim = -1.3,
            L_vplsw = true)
        reecb = Library.reec_b_pf(
            V_dip   = -99,
            V_up    = 99,
            T_rv    = 0.0001, #TODO T_rv=0 Error
            V_ref0  = 1,
            dbd1    = -0.05,
            dbd2    = 0.05,
            K_qv    = 0.0,
            I_qh1   = 1.05,
            I_ql1   = -1.05,
            T_p     = 0.05,
            Q_min   = -0.436,
            Q_max   = 0.436,
            V_min   = 0.9,
            V_max   = 1.1,
            K_qp    = 0,
            K_qi    = 0.1,
            K_vp    = 0,
            K_vi    = 40,
            I_max   = 1.82,
            T_iq    = 0.02,
            T_pord  = 0.02,
            P_min   = 0.0,
            P_max   = 1.0,
            dP_min  = -99,
            dP_max  = 99,
            PqFlag  = false,
            QFlag   = false,
            PfFlag  = false,
            Vflag   = false)
        repca = Library.repc_a_pf(
            K_p        = 18,
            K_i        = 5,
            T_fltr1    = 0.02,
            T_fltr2    = 0.02,
            T_ft       = 0,
            T_fv       = 0.075,
            V_frz      = 0,
            R_c        = 0.0025,
            X_c        = 0.0025,
            K_c        = 0.02,
            e_max      = 0.1,
            e_min      = -0.1,
            dbd_up     = 0.0,
            dbd_dn     = 0.0,
            Q_max      = 0.436,
            Q_min      = -0.436,
            K_pg       = 0.1,
            K_ig       = 0.05,
            T_p        = 0.25,
            fdbd1      = 0,
            fdbd2      = 0,
            femax      = 999,
            femin      = -999,
            P_max      = 999,
            P_min      = -999,
            T_lag      = 0.1,
            D_dn       = 20.0,
            D_up       = 0.0,
            P_plantref = 0.015,
            freq_ref   = 60.0,
            RefFlag    = false,
            VcombFlag  = false,
            freqFlag   = false)
        f = Blocks.Constant(k=60.0)
        Vref = Blocks.Constant(k=1.0)
        Qref = Blocks.Constant(k=-0.056658)
        #Qbranch = Blocks.Constant(k=-0.0567) #TODO nicht konstant
        #Pbranch = Blocks.Constant(k=0.015) #TODO nicht konstant
        #Vreg = Blocks.Constant(k=1.0) #TODO nicht konstant
        Pfaref = Blocks.Constant(k=-1.31199)
        #Vdiff = Blocks.Constant(k=1) #TODO nicht konstant
    end
    @variables begin
        V_t(t), [guess=1, description="Raw terminal voltage"]
        δ_v(t), [guess=0.00045, description="voltage angle"]
        pir(t), [guess=-0.014974521, description="negative terminal current real part"]
        pii(t), [guess=-0.056663541, description="negative terminal current im part"]
        pvr(t), [guess=0.9999999, description=""]
        pvi(t), [guess=0.00044967497, description=""]
        P_gen(t), [guess=0.015, description=""]
        Q_gen(t), [guess=-0.056656801, description=""]
        Vdiff(t), [guess=1.0001042, description=""]
        Vreg_re(t), [guess=1, description=""]
        Vreg_im(t), [guess=0, description=""]
        Qbranch(t), [guess=-0.056656801, description=""]
        Pbranch(t), [guess=0.015, description=""]
    end
    @parameters begin
        S_b=100e6, [description="System Base Power"]
        #M_b=S_b, [description="Machine Base Power"]
        #R_c=0.0025, [description="Line drop compensation resistance"]
        #X_c=0.0025, [description="Line drop compensation reactance"]
    end
    @equations begin
        V_t ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        δ_v ~ atan(terminal.u_i, terminal.u_r)
        regca.Vt_in.u ~ V_t
        reecb.Vt_in.u ~ V_t
        pii ~ - terminal.i_i
        pir ~ - terminal.i_r
        pvr ~ terminal.u_r
        pvi ~ terminal.u_i
        P_gen ~ -(1/CoB) * (pvr*pir + pvi*pii)
        Q_gen ~ -(1/CoB) * (pvi*pir - pvr*pii)
        Pbranch ~ (1/CoB) * (pvr*terminal.i_r + pvi*terminal.i_i)
        Qbranch ~ (1/CoB) * (pvi*terminal.i_r - pvr*terminal.i_i)
        Vreg_re ~ terminal.u_r
        Vreg_im ~ terminal.u_i
        #Vdiff ~ sqrt((terminal.u_r - R_c*terminal.i_r/CoB + X_c*terminal.i_i/CoB)^2 + (terminal.u_i - X_c*terminal.i_r/CoB - R_c*terminal.i_i/CoB)^2)
        #Current injection from converter to terminal (negative for generation)
        [regca.Ipout.u, regca.Iqout.u].~ -1/CoB*[cos(δ_v) sin(δ_v); -sin(δ_v) cos(δ_v)] * [pir, pii]
        connect(repca.freq, f.output)
        connect(repca.V_ref, Vref.output)
        connect(repca.Q_ref, Qref.output)
        repca.Q_branch.u ~ Qbranch
        repca.P_branch.u ~ Pbranch
        repca.V_reg_re.u ~ Vreg_re
        repca.V_reg_im.u ~ Vreg_im
        repca.I_branch.u ~ sqrt(terminal.i_r^2+terminal.i_i^2) #TODO
        reecb.P_e.u ~ P_gen
        reecb.Q_gen.u ~ Q_gen
        reecb.P_inital.u ~ P_gen #TODO
        reecb.Q_initial.u ~ Q_gen #TODO
        connect(reecb.P_faref, Pfaref.output)
        connect(repca.Pref_out, reecb.Pref_in)
        connect(repca.Qext_out, reecb.Qext_in)
        connect(reecb.Iqcmd_out, regca.Iqcmd_in)
        connect(reecb.Ipcmd_out, regca.Ipcmd_in)
    end
end
