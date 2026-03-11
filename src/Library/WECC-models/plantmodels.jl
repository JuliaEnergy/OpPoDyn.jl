# WECC Model for  large scale PV plant taken from OpenIPSL Modelica PSSE (or https://www.esig.energy/wiki-main-page/generic-models-pv-plants/#WECC_Generic_Model_for_Large-scale_PV_Plants)
# default parameter taken from OpenIPSL Examples

@mtkmodel WECC_large_PV begin
    @components begin
        terminal=Terminal()
        regca = regc_a(
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
        reecb = Library.reec_b(
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
        repca = Library.repc_a(
            K_p        = 18,
            K_i        = 5,
            T_fltr     = 0.02,
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
            freqFlag   = false,
            p_0        = 0.015)
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
        Vreg(t), [guess=1, description=""]
        Qbranch(t), [guess=-0.056656801, description=""]
        Pbranch(t), [guess=0.015, description=""]
    end
    @parameters begin
        S_b=100e6, [description="System Base Power"]
        M_b=S_b, [description="Machine Base Power"]
        R_c=0.0025, [description="Line drop compensation resistance"]
        X_c=0.0025, [description="Line drop compensation reactance"]
    end
    begin
        #derived parameters
        CoB = M_b/S_b
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
        Vreg ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        Vdiff ~ sqrt((terminal.u_r - R_c*terminal.i_r/CoB + X_c*terminal.i_i/CoB)^2 + (terminal.u_i - X_c*terminal.i_r/CoB - R_c*terminal.i_i/CoB)^2)
        #Current injection from converter to terminal (negative for generation)
        [regca.Ipout.u, regca.Iqout.u].~ -1/CoB*[cos(δ_v) sin(δ_v); -sin(δ_v) cos(δ_v)] * [pir, pii]
        connect(repca.freq, f.output)
        connect(repca.V_ref, Vref.output)
        connect(repca.Q_ref, Qref.output)
        #connect(repca.Q_branch, Qbranch.output)
        #connect(repca.P_branch, Pbranch.output)
        #connect(repca.V_reg, Vreg.output)
        #connect(repca.V_diff, Vdiff.output)
        repca.Q_branch.u ~ Qbranch
        repca.P_branch.u ~ Pbranch
        repca.V_reg.u ~ Vreg
        repca.V_diff.u ~ Vdiff
        #connect(repca.I_branch, Ibranch.output)
        #connect(reecb.P_e, Pe.output)
        #connect(reecb.Q_gen, Qgen.output)
        reecb.P_e.u ~ P_gen
        reecb.Q_gen.u ~ Q_gen
        connect(reecb.P_faref, Pfaref.output)
        #connect(regca.Q_gen0, Qgen0.output)
        connect(repca.Pref_out, reecb.Pref_in)
        connect(repca.Qext_out, reecb.Qext_in)
        connect(reecb.Iqcmd_out, regca.Iqcmd_in)
        connect(reecb.Ipcmd_out, regca.Ipcmd_in)
    end
end


@mtkmodel WECC_BESS begin
    @components begin
        terminal=Terminal()
        regca = regc_a(
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
        reecc = Library.reec_c(
            V_dip   = -99,
            V_up    = 99,
            T_rv    = 0.01,
            V_ref0  = 1,
            dbd1    = -0.05,
            dbd2    = 0.05,
            K_qv    = 15,
            I_qh1   = 0.75,
            I_ql1   = -0.75,
            T_p     = 0.05,
            Q_min   = -0.75,
            Q_max   = 0.75,
            V_min   = 0.9,
            V_max   = 1.1,
            K_qp    = 0,
            K_qi    = 1,
            K_vp    = 0,
            K_vi    = 1,
            I_max   = 1.11,
            T_iq    = 0.017,
            T_pord  = 0.017,
            P_min   = -0.667,
            P_max   = 1.0,
            dP_min  = -99,
            dP_max  = 99,
            soc_ini = 0.5,
            T_char = 999,
            SOCmin = 0.2,
            SOCmax = 0.8,
            PfFlag  = false,
            Vflag   = false,
            QFlag   = false,
            PqFlag  = false)
        repca = Library.repc_a(
            K_p        = 18,
            K_i        = 5,
            T_fltr     = 0.02,
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
            P_plantref = 0.015, #from powerflow
            freq_ref   = 50.0,
            p_0        = 0.015,
            RefFlag    = false,
            VcombFlag  = false,
            freqFlag   = false)
        f = Blocks.Constant(k=50.0)
        Vref = Blocks.Constant(k=1.0)
        Qref = Blocks.Constant(k=-0.056658)
        #Qbranch = Blocks.Constant(k=-0.0567)
        #Pbranch = Blocks.Constant(k=0.015)
        #Vreg = Blocks.Constant(k=1.0)
        #Ibranch = Blocks.Constant(k=1.0)
        #Pe = Blocks.Constant(k=1.0)
        Pfaref = Blocks.Constant(k=-1.31199)
        #Qgen = Blocks.Constant(k=0.0)
        #Qgen0 = Blocks.Constant(k=0.0)
        #Vdiff = Blocks.Constant(k=1)
        P_aux = Blocks.Constant(k=0)
        #PELEC = Blocks.Constant(k=0.015)
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
        PELEC(t), [guess=0.015, description=""]
        Vdiff(t), [guess=1.0001042, description=""]
        Vreg(t), [guess=1, description=""]
        Qbranch(t), [guess=-0.056656801, description=""]
        Pbranch(t), [guess=0.015, description=""]
    end
    @parameters begin
        S_b=100e6, [description="System Base Power"]
        M_b=S_b, [description="Machine Base Power"]
        R_c=0.0025, [description="Line drop compensation resistance"]
        X_c=0.0025, [description="Line drop compensation reactance"]
    end
    begin
        #derived parameters
        CoB = M_b/S_b
    end
    @equations begin
        V_t ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        δ_v ~ atan(terminal.u_i, terminal.u_r)
        regca.Vt_in.u ~ V_t
        reecc.Vt_in.u ~ V_t
        pii ~ - terminal.i_i
        pir ~ - terminal.i_r
        pvr ~ terminal.u_r
        pvi ~ terminal.u_i
        P_gen ~ -(1/CoB) * (pvr*pir + pvi*pii)
        Q_gen ~ -(1/CoB) * (pvi*pir - pvr*pii)
        Pbranch ~ (1/CoB) * (pvr*terminal.i_r + pvi*terminal.i_i)
        Qbranch ~ (1/CoB) * (pvi*terminal.i_r - pvr*terminal.i_i)
        Vreg ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        Vdiff ~ sqrt((terminal.u_r - R_c*terminal.i_r/CoB + X_c*terminal.i_i/CoB)^2 + (terminal.u_i - X_c*terminal.i_r/CoB - R_c*terminal.i_i/CoB)^2)
        [regca.Ipout.u, regca.Iqout.u].~ -1/CoB*[cos(δ_v) sin(δ_v); -sin(δ_v) cos(δ_v)] * [pir, pii]
        connect(repca.freq, f.output)
        connect(repca.V_ref, Vref.output)
        connect(repca.Q_ref, Qref.output)
        #connect(repca.V_diff, Vdiff.output)
        #connect(repca.Q_branch, Qbranch.output)
        #connect(repca.P_branch, Pbranch.output)
        #connect(repca.V_reg, Vreg.output)
        repca.Q_branch.u ~ Qbranch
        repca.P_branch.u ~ Pbranch
        repca.V_reg.u ~ Vreg
        repca.V_diff.u ~ Vdiff
        #connect(repca.I_branch, Ibranch.output)
        #connect(reecc.P_e, Pe.output)
        connect(reecc.P_faref, Pfaref.output)
        #connect(reecc.Q_gen, Qgen.output)
        reecc.P_e.u ~ P_gen
        reecc.Q_gen.u ~ Q_gen
        #connect(reecc.PELEC, PELEC.output)
        reecc.PELEC.u ~ P_gen
        connect(reecc.P_aux, P_aux.output)
        #connect(regca.Q_gen0, Qgen0.output)
        connect(repca.Pref_out, reecc.Pref_in)
        connect(repca.Qext_out, reecc.Qext_in)
        connect(reecc.Iqcmd_out, regca.Iqcmd_in)
        connect(reecc.Ipcmd_out, regca.Ipcmd_in)
    end
end


@mtkmodel WECC_WT_4B begin
    @components begin
        terminal=Terminal()
        regca = regc_a(
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
        reeca = Library.reec_a(
            V_0     = 1,
            V_dip   = -99, #0.85,
            V_up    = 99, #1.2,
            T_rv    = 0.01,
            V_ref0  = 1,
            dbd1    = -0.1,
            dbd2    = 0.1,
            K_qv    = 0.0,
            I_qh1   = 1,
            I_ql1   = -1,
            T_p     = 0.3,
            Q_min   = -0.5,
            Q_max   = 0.5,
            V_min   = 0.9,
            V_max   = 1.1,
            K_qp    = 0.1,
            K_qi    = 0.1,
            V_bias  = 0,
            K_vp    = 1.6,
            K_vi    = 1,
            I_max   = 1.2,
            T_iq    = 0.01,
            T_pord  = 0.3,
            P_min   = 0.0,
            P_max   = 1.0,
            dP_min  = -99,
            dP_max  = 99,
            PfFlag  = false,
            Vflag   = false,
            QFlag   = false,
            PqFlag  = false)
        repca = Library.repc_a(
            K_p        = 18,
            K_i        = 5,
            T_fltr     = 0.02,
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
            P_plantref = 0.015, #from powerflow
            freq_ref   = 50.0,
            p_0        = 0.015,
            RefFlag    = false,
            VcombFlag  = false,
            freqFlag   = false)
        drive_train = Library.WTDTA1()
        f = Blocks.Constant(k=50.0)
        Vref = Blocks.Constant(k=1.0)
        Qref = Blocks.Constant(k=-0.056658)
        #Qbranch = Blocks.Constant(k=-0.0567)
        #Pbranch = Blocks.Constant(k=0.015)
        #Vreg = Blocks.Constant(k=1.0)
        #Ibranch = Blocks.Constant(k=1.0)
        #Pe = Blocks.Constant(k=1.0)
        Pfaref = Blocks.Constant(k=-1.31199)
        #Qgen = Blocks.Constant(k=0.0)
        #Qgen0 = Blocks.Constant(k=0.0)
        #Vdiff = Blocks.Constant(k=1)
        #WG = Blocks.Constant(k=1)
        W0 = Blocks.Constant(k=0)
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
        Vreg(t), [guess=1, description=""]
        Qbranch(t), [guess=-0.056656801, description=""]
        Pbranch(t), [guess=0.015, description=""]
    end
    @parameters begin
        S_b=100e6, [description="System Base Power in VA"]
        M_b=S_b, [description="Machine Base Power in VA"]
        P_0=1.5e6, [description="Initial activa power in W"]
        R_c=0.0025, [description="Line drop compensation resistance"]
        X_c=0.0025, [description="Line drop compensation reactance"]
    end
    begin
        #derived parameters
        CoB = M_b/S_b
        p0 = P_0/M_b
    end
    @equations begin
        V_t ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        δ_v ~ atan(terminal.u_i, terminal.u_r)
        pii ~ - terminal.i_i
        pir ~ - terminal.i_r
        pvr ~ terminal.u_r
        pvi ~ terminal.u_i
        [regca.Ipout.u, regca.Iqout.u].~ -1/CoB*[cos(δ_v) sin(δ_v); -sin(δ_v) cos(δ_v)] * [pir, pii]
        P_gen ~ -(1/CoB) * (pvr*pir + pvi*pii)
        Q_gen ~ -(1/CoB) * (pvi*pir - pvr*pii)
        Pbranch ~ (1/CoB) * (pvr*terminal.i_r + pvi*terminal.i_i)
        Qbranch ~ (1/CoB) * (pvi*terminal.i_r - pvr*terminal.i_i)
        Vreg ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        Vdiff ~ sqrt((terminal.u_r - R_c*terminal.i_r/CoB + X_c*terminal.i_i/CoB)^2 + (terminal.u_i - X_c*terminal.i_r/CoB - R_c*terminal.i_i/CoB)^2)
        regca.Vt_in.u ~ V_t
        reeca.Vt_in.u ~ V_t
        drive_train.P_e.u ~ P_gen
        drive_train.P_m.u ~ p0
        #connect(regca.Q_gen0, Qgen0.output)
        connect(repca.freq, f.output)
        connect(repca.V_ref, Vref.output)
        connect(repca.Q_ref, Qref.output)
        #connect(repca.Q_branch, Qbranch.output)
        #connect(repca.P_branch, Pbranch.output)
        #connect(repca.V_reg, Vreg.output)
        #connect(repca.V_diff, Vdiff.output)
        repca.Q_branch.u ~ Qbranch
        repca.P_branch.u ~ Pbranch
        repca.V_reg.u ~ Vreg
        repca.V_diff.u ~ Vdiff
        #connect(repca.I_branch, Ibranch.output)
        #connect(reeca.P_e, Pe.output)
        connect(reeca.P_faref, Pfaref.output)
        #connect(reeca.Q_gen, Qgen.output)
        reeca.P_e.u ~ P_gen
        reeca.Q_gen.u ~ Q_gen
        connect(drive_train.W_0, W0.output)
        connect(drive_train.W_gout, reeca.Wg)
        #connect(reeca.Wg, WG.output)
        connect(repca.Pref_out, reeca.Pref_in)
        connect(repca.Qext_out, reeca.Qext_in)
        connect(reeca.Iqcmd_out, regca.Iqcmd_in)
        connect(reeca.Ipcmd_out, regca.Ipcmd_in)
    end
end


@mtkmodel WTDTA1 begin
    @components begin
        P_m = RealInput(guess=0.015)
        P_e = RealInput(guess=0.015)
        W_0 = RealInput(guess=0)
        W_tout = RealOutput(guess=0)
        W_gout = RealOutput(guess=1)
    end
    @parameters begin
        D_shaft=1, [description="Shaft damping factor"]
        H=5.3, [description="Total inertia constant"]
        H_tfrac=0.92, [description="Turbine inertia fraction, H_t/H"]
        K_damp=0, [description=""]
        fn=50, [description="System nominal frequency in Hz"]
        freq1=2.132, [description="First shaft torsional resonancy frequency"]
    end
    @variables begin
        w_tin(t), [guess=0, description=""]
        w_add(t), [guess=0.015, description=""]
        w_tsum(t), [guess=1, description=""]
        w_t(t), [guess=0, description=""]
        w_gin(t), [guess=0, description=""]
        Δw_g(t), [guess=0, description=""]
        w_g(t), [guess=1, description=""]
        w_gint(t), [guess=0, description=""]
        Δw(t), [guess=0, description=""]
        #w_add3(t), [description=""] #not connected
        #Δw_g3(t), [description=""]
    end
    begin
        H_t = H * H_tfrac
        H_g = H - H_t
        w0 = 2 * π * fn
        K_shaft= 2 * H_t * H_g * (2 * π * freq1)^2/(H * w0)
    end
    @equations begin
        w_tin ~ P_m.u/w_tsum - w_add - D_shaft * Δw
        2 * H_t * Dt(w_t) ~ w_tin
        w_tsum ~ w_t + 1
        w_gin ~ D_shaft * Δw + w_add - P_e.u/w_g
        Δw_g ~ w_gin - K_damp * w_gint
        2 * H_g * Dt(w_gint) ~ Δw_g
        w_g ~ w_gint + 1
        Δw ~ w_t - w_gint
        Dt(w_add) ~ w0 * K_shaft * Δw
        #Δw_g3 ~ w_gint - W_0.u
        #Dt(w_gadd3) ~ w0 * Δw_g3
        W_gout.u ~ w_g
        W_tout.u ~ w_t
    end
end


@mtkmodel WECC_large_PV_pf begin
    @structural_parameters begin
        f_input = false
        Vref_input = false
        Qref_input = false
        Qinit_input = false
        Pinit_input = false
        P_plantref_input = false
    end
    @components begin
        terminal=Terminal()
        regca = regc_a_pf(
            I_qrmax = 999,
            I_qrmin = -999,
            T_gp    = 0.02,
            T_gq    = 0.02,
            T_fltr  = 0.02,
            Brkpt   = 0.9,
            Zerox   = 0.4,
            L_vpl1  = 1.22,
            rrpwr   = 10,
            V_0lim = 1.2,
            K_hv = 0.7,
            lvpnt0 = 0.4,
            lvpnt1 = 0.8,
            I_0lim = -1.3,
            L_vplsw = false)
        reecb = Library.reec_b_pf(
            V_dip   = 0,
            V_up    = 999,
            T_rv    = 0.02,
            V_ref0  = 0,
            dbd1    = -0.02,
            dbd2    = 0.02,
            K_qv    = 0.0,
            I_qh1   = 1.05,
            I_ql1   = -1.05,
            T_p     = 0.02,
            Q_min   = -0.4,
            Q_max   = 0.4,
            V_min   = 0.9,
            V_max   = 1.1,
            K_qp    = 1,
            K_qi    = 10,
            K_vp    = 0,
            K_vi    = 40,
            I_max   = 1.25,
            T_iq    = 0.02,
            T_pord  = 0.02,
            P_min   = 0.0,
            P_max   = 1.0,
            dP_min  = -999,
            dP_max  = 999,
            PqFlag  = false,
            QFlag   = false,
            PfFlag  = false,
            Vflag   = true)
        repca = Library.repc_a_pf(
            K_p        = 1,
            K_i        = 10,
            T_fltr1     = 0.02,
            T_fltr2     = 0.02,
            T_ft       = 0,
            T_fv       = 0.1,
            V_frz      = 0,
            R_c        = 0,
            X_c        = 0,
            K_c        = 0.02,
            e_max      = 0.1,
            e_min      = -0.1,
            dbd_up     = 0.0,
            dbd_dn     = 0.0,
            Q_max      = 0.4,
            Q_min      = -0.4,
            K_pg       = 1,
            K_ig       = 10,
            T_p        = 0.2,
            fdbd1      = 0,
            fdbd2      = 0,
            femax      = 999,
            femin      = -999,
            P_max      = 1,
            P_min      = 0,
            T_lag      = 0.1,
            D_dn       = 20.0,
            D_up       = 10.0,
            freq_ref   = 1, #in p.u.
            RefFlag    = false,
            VcombFlag  = true,
            freqFlag   = false)
        if f_input
            f_in = RealInput()
        end
        if Vref_input
            Vref_in = RealInput()
        end
        if Qref_input
            Qref_in = RealInput()
        end
        if Qinit_input
            Qinit_in = RealInput()
        end
        if Pinit_input
            Pinit_in = RealInput()
        end
        if P_plantref_input
            P_plantref_in = RealInput()
        end
    end
    @variables begin
        V_t(t), [guess=1.001047, description="Raw terminal voltage"]
        δ_v(t), [guess=0.0023076, description="voltage angle in rad"]
        pir(t), [guess=0.799189, description="negative terminal current real part"]
        pii(t), [guess=0.301802, description="negative terminal current im part"]
        pvr(t), [guess=1.001044, description=""]
        pvi(t), [guess=0.00231, description=""]
        #P_gen(t), [guess=0.800721, description=""]
        #Q_gen(t), [guess=-0.30027, description=""]
        #Vdiff(t), [guess=, description=""]
        #Vreg_re(t), [guess=1.001044, description=""]
        #Vreg_im(t), [guess=0.00231, description=""]
        Q_measure(t), [guess=-0.30027, description=""]
        P_measure(t), [guess=0.800721, description=""]
        I_measure(t), [guess=0.854276, description=""]
        P_plantref(t), [guess=0.800721, description=""]
        f(t), [guess=1, description=""]
        V_ref(t), [guess=1.001047, description=""]
        Q_ref(t), [guess=-0.30027, description=""]
        Qinit(t), [guess=-0.30027, description=""]
        Pinit(t), [guess=0.800721, description=""]
    end
    @parameters begin
        S_b=100e6, [description="System Base Power"]
        #M_b=S_b, [description="Machine Base Power"]
        if !f_input
            f_set=1, [guess=1, description="frequency setpoint [pu]"]
        end
        if !Vref_input
            Vref_set=1.001046, [guess=1.001046, description="[pu]"]
        end
        if !Qref_input
            Qref_set=-0.30027, [guess=-0.30027, description="[pu]"]
        end
        if !Qinit_input
            Qinit_set=-0.30027, [guess=-0.30027, description="[pu]"]
        end
        if !Pinit_input
            Pinit_set=0.800721, [guess=0.800721, description="[pu]"]
        end
        if !P_plantref_input
            P_plantref_set=0.800721, [guess=0.800721, description="[pu]"]
        end
    end
    @equations begin
        f ~ f_input ? f_in.u : f_set
        V_ref ~ Vref_input ? Vref_in.u : Vref_set
        Q_ref ~ Qref_input ? Qref_in.u : Qref_set
        Qinit ~ Qinit_input ? Qinit_in.u : Qinit_set
        Pinit ~ Pinit_input ? Pinit_in.u : Pinit_set
        P_plantref ~ P_plantref_input ? P_plantref_in.u : P_plantref_set
        #calculate inputs
        V_t ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        δ_v ~ atan(terminal.u_i, terminal.u_r)
        pii ~ - terminal.i_i
        pir ~ - terminal.i_r
        pvr ~ terminal.u_r
        pvi ~ terminal.u_i
        #P_gen ~ - (pvr*pir + pvi*pii)
        #Q_gen ~ - (pvi*pir - pvr*pii)
        I_measure ~ sqrt(pii^2 + pir^2)
        P_measure ~  (pvr*terminal.i_r + pvi*terminal.i_i)
        #P_plantref ~ (pvr*terminal.i_r + pvi*terminal.i_i) #TODO Wert in PF
        Q_measure ~  (pvi*terminal.i_r - pvr*terminal.i_i)
        #Vreg_re ~ terminal.u_r
        #Vreg_im ~ terminal.u_i
        #Current injection from converter to terminal (negative for generation)
        [regca.Ipout.u, regca.Iqout.u].~ -[cos(δ_v) sin(δ_v); -sin(δ_v) cos(δ_v)] * [pir, pii] #TODO entspricht das PF Berechnung?
        #[regca.Ipout.u, regca.Iqout.u].~ -[cos(0) sin(0); -sin(0) cos(0)] * [pir, pii] #TODO entspricht das PF Berechnung?
        #[regca.Ipout.u, regca.Iqout.u] .~ [pir, pii]
        #regca.Ipout.u ~ pir
        #regca.Iqout.u ~ pii

        #connect components
        repca.freq.u ~ f
        repca.V_ref.u ~ V_ref
        repca.Q_ref.u ~ Q_ref
        repca.P_plantref.u ~ P_plantref 
        repca.Q_branch.u ~ Q_measure
        repca.P_branch.u ~ P_measure
        repca.I_branch.u ~ I_measure
        repca.V_reg_re.u ~ pvr
        repca.V_reg_im.u ~ pvi

        reecb.Vt_in.u ~ V_t
        reecb.P_e.u ~ P_measure
        reecb.Q_gen.u ~ Q_measure
        reecb.P_initial.u ~ Pinit
        reecb.Q_initial.u ~ Qinit

        connect(repca.Pref_out, reecb.Pref_in)
        connect(repca.Qext_out, reecb.Qext_in)
        connect(reecb.Iqcmd_out, regca.Iqcmd_in)
        connect(reecb.Ipcmd_out, regca.Ipcmd_in)
        connect(regca.V_tfilt, reecb.Vt_filt)
    end
end

