# WECC Model for  large scale PV plant taken from OpenIPSL Modelica PSSE (or https://www.esig.energy/wiki-main-page/generic-models-pv-plants/#WECC_Generic_Model_for_Large-scale_PV_Plants)
# default parameter taken from OpenIPSL Examples

@mtkmodel WECC_large_PV begin
    @components begin
        terminal=Terminal()
        converter_interface = regc_a(
            I_qrmax = 9999,
            I_qrmin = -9999,
            T_g     = 0.02,
            T_fltr  = 0.02,
            Brkpt   = 0.9,
            Zerox   = 0.5,
            L_vpl1  = 1.22,
            L_vplsw = 1,
            rrpwr   = 10,
            V_0lim = 1.2,
            K_hv = 0.7,
            lvpnt0 = 0.4,
            lvpnt1 = 0.8,
            I_0lim = -1.3)
        electrical_control = Library.reec_b(
            V_dip   = -99, #0.85
            V_up    = 99, #1.2
            T_rv    = 0.01,
            V_ref0  = 1,
            dbd1    = -0.1,
            dbd2    = 0.1,
            K_qv    = 0.0,
            I_qh1   = 1.05,
            I_ql1   = -1.05,
            T_p     = 0.051,
            PfFlag  = 0,
            Q_min   = -0.436,
            Q_max   = 0.436,
            V_min   = 0.9,
            V_max   = 1.1,
            K_qp    = 0,
            K_qi    = 0.1,
            Vflag   = 0,
            K_vp    = 0,
            K_vi    = 40,
            QFlag   = 0,
            I_max   = 1.82,
            PqFlag  = 0,
            T_iq    = 0.02,
            T_pord  = 0.02,
            P_min   = 0.0,
            P_max   = 1.0,
            dP_min  = -99,
            dP_max  = 99)
        plant_control = Library.repc_a(
            K_p        = 18,
            K_i        = 5,
            T_fltr     = 0.02,
            T_ft       = 0,
            T_fv       = 0.075,
            RefFlag    = 0,
            V_frz      = 0,
            R_c        = 0.0025,
            X_c        = 0.0025,
            K_c        = 0.02,
            VcombFlag  = 1,
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
            freq_ref   = 60.0,
            freqFlag   = 1,
            p_0        = 0.015)
        f = Blocks.Constant(k=60.0)
        Vref = Blocks.Constant(k=1.0)
        Qref = Blocks.Constant(k=-0.0567)
        Qbranch = Blocks.Constant(k=-0.0567)
        Pbranch = Blocks.Constant(k=0.015)
        Vreg = Blocks.Constant(k=1.0)
        #Ibranch = Blocks.Constant(k=1.0)
        #Pe = Blocks.Constant(k=1.0)
        Pfaref = Blocks.Constant(k=-1.31199)
        #Qgen = Blocks.Constant(k=0.0)
        #Qgen0 = Blocks.Constant(k=0.0)
        Vdiff = Blocks.Constant(k=1)
    end
    @variables begin
        V_t(t), [guess=1, description="Raw terminal voltage"]
        δ_v(t), [guess=0.00045, description="voltage angle"]
        pir(t), [guess=-0.015, description="negative terminal current real part"]
        pii(t), [guess=-0.0567, description="negative terminal current im part"]
        pvr(t), [guess=1, description=""]
        pvi(t), [guess=0.00045, description=""]
        P_gen(t), [guess=0.015, description=""]
        Q_gen(t), [guess=-0.0567, description=""]
    end
    @parameters begin
        S_b=100e6, [description="System Base Power"]
        M_b=S_b, [description="Machine Base Power"]
    end
    begin
        #derived parameters
        CoB = M_b/S_b
    end
    @equations begin
        V_t ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        δ_v ~ atan(terminal.u_i, terminal.u_r)
        converter_interface.Vt_in.u ~ V_t
        electrical_control.Vt_in.u ~ V_t
        pii ~ - terminal.i_i
        pir ~ - terminal.i_r
        pvr ~ terminal.u_r
        pvi ~ terminal.u_i
        P_gen ~ -(1/CoB) * (pvr*pir + pvi*pii)
        Q_gen ~ -(1/CoB) * (pvi*pir - pvr*pii)
        #Current injection from converter to terminal (negative for generation)
        [converter_interface.Ipout.u, converter_interface.Iqout.u].~ -1/CoB*[cos(δ_v) sin(δ_v); -sin(δ_v) cos(δ_v)] * [pir, pii]
        connect(plant_control.V_diff, Vdiff.output)
        connect(plant_control.freq, f.output)
        connect(plant_control.V_ref, Vref.output)
        connect(plant_control.Q_ref, Qref.output)
        connect(plant_control.Q_branch, Qbranch.output)
        connect(plant_control.P_branch, Pbranch.output)
        connect(plant_control.V_reg, Vreg.output)
        #connect(plant_control.I_branch, Ibranch.output)
        #connect(electrical_control.P_e, Pe.output)
        #connect(electrical_control.Q_gen, Qgen.output)
        electrical_control.P_e.u ~ P_gen
        electrical_control.Q_gen.u ~ Q_gen
        connect(electrical_control.P_faref, Pfaref.output)
        #connect(converter_interface.Q_gen0, Qgen0.output)
        connect(plant_control.Pref_out, electrical_control.Pref_in)
        connect(plant_control.Qext_out, electrical_control.Qext_in)
        connect(electrical_control.Iqcmd_out, converter_interface.Iqcmd_in)
        connect(electrical_control.Ipcmd_out, converter_interface.Ipcmd_in)
    end
end


@mtkmodel WECC_BESS begin
    @components begin
        terminal=Terminal()
        converter_interface = regc_a(
            I_qrmax = 9999,
            I_qrmin = -9999,
            T_g     = 0.02,
            T_fltr  = 0.02,
            Brkpt   = 0.9,
            Zerox   = 0.5,
            L_vpl1  = 1.22,
            L_vplsw = 0,
            rrpwr   = 10,
            V_0lim = 1.2,
            K_hv = 0.7,
            lvpnt0 = 0.4,
            lvpnt1 = 0.8,
            I_0lim = -1.3)
        electrical_control = Library.reec_c(
            V_dip   = -99, #0.85,
            V_up    = 99, #1.2,
            T_rv    = 0.01,
            V_ref0  = 1,
            dbd1    = 0,
            dbd2    = 0,
            K_qv    = 0.0,
            I_qh1   = 1,
            I_ql1   = -1,
            T_p     = 0.01,
            PfFlag  = 0,
            Q_min   = -1,
            Q_max   = 1,
            V_min   = 0.9,
            V_max   = 1.1,
            K_qp    = 0,
            K_qi    = 0.1,
            Vflag   = 0,
            K_vp    = 0,
            K_vi    = 40,
            QFlag   = 0,
            I_max   = 1.11,
            PqFlag  = 0,
            T_iq    = 0.01,
            T_pord  = 0.02,
            P_min   = 0.0,
            P_max   = 1.0,
            dP_min  = -99,
            dP_max  = 99,
            soc_ini = 0.5,
            T_char = 999,
            SOCmin = 0.2,
            SOCmax = 0.8)
        plant_control = Library.repc_a(
            K_p        = 18,
            K_i        = 5,
            T_fltr     = 0.02,
            T_ft       = 0,
            T_fv       = 0.075,
            RefFlag    = 0,
            V_frz      = 0,
            R_c        = 0.0025,
            X_c        = 0.0025,
            K_c        = 0.02,
            VcombFlag  = 1,
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
            freqFlag   = 1,
            p_0        = 0.015)
        f = Blocks.Constant(k=50.0)
        Vref = Blocks.Constant(k=1.0)
        Qref = Blocks.Constant(k=-0.0567)
        Qbranch = Blocks.Constant(k=-0.0567)
        Pbranch = Blocks.Constant(k=0.015)
        Vreg = Blocks.Constant(k=1.0)
        #Ibranch = Blocks.Constant(k=1.0)
        #Pe = Blocks.Constant(k=1.0)
        Pfaref = Blocks.Constant(k=-1.311)
        #Qgen = Blocks.Constant(k=0.0)
        #Qgen0 = Blocks.Constant(k=0.0)
        Vdiff = Blocks.Constant(k=1)
        P_aux = Blocks.Constant(k=0)
        PELEC = Blocks.Constant(k=0.015)
    end
    @variables begin
        V_t(t), [guess=1, description="Raw terminal voltage"]
        δ_v(t), [guess=0.00045, description="voltage angle in rad"]
        pir(t), [guess=-0.015, description="negative terminal current real part"]
        pii(t), [guess=-0.0567, description="negative terminal current im part"]
        pvr(t), [guess=1, description=""]
        pvi(t), [guess=0.00045, description=""]
        P_gen(t), [guess=0.015, description=""]
        Q_gen(t), [guess=-0.0567, description=""]
    end
    @parameters begin
        S_b=100e6, [description="System Base Power"]
        M_b=S_b, [description="Machine Base Power"]
    end
    begin
        #derived parameters
        CoB = M_b/S_b
    end
    @equations begin
        V_t ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        δ_v ~ atan(terminal.u_i, terminal.u_r)
        converter_interface.Vt_in.u ~ V_t
        electrical_control.Vt_in.u ~ V_t
        pii ~ - terminal.i_i
        pir ~ - terminal.i_r
        pvr ~ terminal.u_r
        pvi ~ terminal.u_i
        P_gen ~ -(1/CoB) * (pvr*pir + pvi*pii)
        Q_gen ~ -(1/CoB) * (pvi*pir - pvr*pii)
        [converter_interface.Ipout.u, converter_interface.Iqout.u].~ -1/CoB*[cos(δ_v) sin(δ_v); -sin(δ_v) cos(δ_v)] * [pir, pii]
        connect(plant_control.V_diff, Vdiff.output)
        connect(plant_control.freq, f.output)
        connect(plant_control.V_ref, Vref.output)
        connect(plant_control.Q_ref, Qref.output)
        connect(plant_control.Q_branch, Qbranch.output)
        connect(plant_control.P_branch, Pbranch.output)
        connect(plant_control.V_reg, Vreg.output)
        #connect(plant_control.I_branch, Ibranch.output)
        #connect(electrical_control.P_e, Pe.output)
        connect(electrical_control.P_faref, Pfaref.output)
        #connect(electrical_control.Q_gen, Qgen.output)
        electrical_control.P_e.u ~ P_gen
        electrical_control.Q_gen.u ~ Q_gen
        connect(electrical_control.PELEC, PELEC.output)
        connect(electrical_control.P_aux, P_aux.output)
        #connect(converter_interface.Q_gen0, Qgen0.output)
        connect(plant_control.Pref_out, electrical_control.Pref_in)
        connect(plant_control.Qext_out, electrical_control.Qext_in)
        connect(electrical_control.Iqcmd_out, converter_interface.Iqcmd_in)
        connect(electrical_control.Ipcmd_out, converter_interface.Ipcmd_in)
    end
end


@mtkmodel WECC_WT_4B begin
    @components begin
        terminal=Terminal()
        converter_interface = regc_a(
            I_qrmax = 9999,
            I_qrmin = -9999,
            T_g     = 0.02,
            T_fltr  = 0.02,
            Brkpt   = 0.9,
            Zerox   = 0.5,
            L_vpl1  = 1.22,
            L_vplsw = 0,
            rrpwr   = 10,
            V_0lim = 1.2,
            K_hv = 0.7,
            lvpnt0 = 0.4,
            lvpnt1 = 0.8,
            I_0lim = -1.3)
        electrical_control = Library.reec_a(
            V_dip   = -99, #0.85,
            V_up    = 99, #1.2,
            T_rv    = 0.01,
            V_ref0  = 1,
            dbd1    = -0.05,
            dbd2    = 0.05,
            K_qv    = 0.0,
            I_qh1   = 1.05,
            I_ql1   = -1.05,
            T_p     = 0.05,
            PfFlag  = 0,
            Q_min   = -0.4,
            Q_max   = 0.4,
            V_min   = 0.9,
            V_max   = 1.1,
            K_qp    = 0,
            K_qi    = 0.1,
            Vflag   = 0,
            V_bias  = 0,
            K_vp    = 1.6,
            K_vi    = 1,
            QFlag   = 0,
            I_max   = 1.2,
            PqFlag  = 0,
            T_iq    = 0.02,
            T_pord  = 0.04,
            P_min   = 0.0,
            P_max   = 1.0,
            dP_min  = -99,
            dP_max  = 99)
        plant_control = Library.repc_a(
            K_p        = 18,
            K_i        = 5,
            T_fltr     = 0.02,
            T_ft       = 0,
            T_fv       = 0.075,
            RefFlag    = 0,
            V_frz      = 0,
            R_c        = 0.0025,
            X_c        = 0.0025,
            K_c        = 0.02,
            VcombFlag  = 1,
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
            freqFlag   = 1,
            p_0        = 0.015)
        drive_train = Library.WTDTA1()
        f = Blocks.Constant(k=50.0)
        Vref = Blocks.Constant(k=1.0)
        Qref = Blocks.Constant(k=-0.0567)
        Qbranch = Blocks.Constant(k=-0.0567)
        Pbranch = Blocks.Constant(k=0.015)
        Vreg = Blocks.Constant(k=1.0)
        #Ibranch = Blocks.Constant(k=1.0)
        #Pe = Blocks.Constant(k=1.0)
        Pfaref = Blocks.Constant(k=-1.311)
        #Qgen = Blocks.Constant(k=0.0)
        #Qgen0 = Blocks.Constant(k=0.0)
        Vdiff = Blocks.Constant(k=1)
        #WG = Blocks.Constant(k=1)
        W0 = Blocks.Constant(k=0)
    end
    @variables begin
        V_t(t), [guess=1, description="Raw terminal voltage"]
        δ_v(t), [guess=0.00045, description="voltage angle"]
        pir(t), [guess=-0.015, description="negative terminal current real part"]
        pii(t), [guess=-0.0567, description="negative terminal current im part"]
        pvr(t), [guess=1, description=""]
        pvi(t), [guess=0.00045, description=""]
        P_gen(t), [guess=0.015, description=""]
        Q_gen(t), [guess=-0.0567, description=""]
    end
    @parameters begin
        S_b=100e6, [description="System Base Power in VA"]
        M_b=S_b, [description="Machine Base Power in VA"]
        P_0=1.5e6, [description="Initial activa power in W"]
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
        [converter_interface.Ipout.u, converter_interface.Iqout.u].~ -1/CoB*[cos(δ_v) sin(δ_v); -sin(δ_v) cos(δ_v)] * [pir, pii]
        P_gen ~ -(1/CoB) * (pvr*pir + pvi*pii)
        Q_gen ~ -(1/CoB) * (pvi*pir - pvr*pii)
        converter_interface.Vt_in.u ~ V_t
        electrical_control.Vt_in.u ~ V_t
        drive_train.P_e.u ~ P_gen
        drive_train.P_m.u ~ p0
        #connect(converter_interface.Q_gen0, Qgen0.output)
        connect(plant_control.V_diff, Vdiff.output)
        connect(plant_control.freq, f.output)
        connect(plant_control.V_ref, Vref.output)
        connect(plant_control.Q_ref, Qref.output)
        connect(plant_control.Q_branch, Qbranch.output)
        connect(plant_control.P_branch, Pbranch.output)
        connect(plant_control.V_reg, Vreg.output)
        #connect(plant_control.I_branch, Ibranch.output)
        #connect(electrical_control.P_e, Pe.output)
        connect(electrical_control.P_faref, Pfaref.output)
        #connect(electrical_control.Q_gen, Qgen.output)
        electrical_control.P_e.u ~ P_gen
        electrical_control.Q_gen.u ~ Q_gen
        connect(drive_train.W_0, W0.output)
        connect(drive_train.W_gout, electrical_control.Wg)
        #connect(electrical_control.Wg, WG.output)
        connect(plant_control.Pref_out, electrical_control.Pref_in)
        connect(plant_control.Qext_out, electrical_control.Qext_in)
        connect(electrical_control.Iqcmd_out, converter_interface.Iqcmd_in)
        connect(electrical_control.Ipcmd_out, converter_interface.Ipcmd_in)
    end
end


@mtkmodel WTDTA1 begin
    @components begin
        P_m = RealInput(guess=0.015)
        P_e = RealInput(guess=0.015)
        W_0 = RealInput(guess=0)
        W_tout = RealOutput()
        W_gout = RealOutput()
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