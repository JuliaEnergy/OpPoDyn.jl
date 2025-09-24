# WECC Model for  large scale PV plant taken from OpenIPSL Modelica PSSE (or https://www.esig.energy/wiki-main-page/generic-models-pv-plants/#WECC_Generic_Model_for_Large-scale_PV_Plants)
# default parameter taken from OpenIPSL, flags=0

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
            L_vplsw = 0,
            rrpwr   = 10,
            V_0lim = 1.2,
            K_hv = 0.7,
            lvpnt0 = 0.4,
            lvpnt1 = 0.8,
            I_0lim = -1.3)
        electrical_control = Library.reec_b(
            V_dip   = 0.85,
            V_up    = 1.2,
            T_rv    = 0.01,
            V_ref0  = 0.95,
            dbd1    = -0.1,
            dbd2    = 0.1,
            K_qv    = 0.0,
            I_qhl   = 1.05,
            I_qll   = -1.05,
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
            VcombFlag  = 0,
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
            P_plantref = 1.0, #from powerflow
            freq_ref   = 1.0,
            freqFlag   = 0)
        f = Blocks.Constant(k=1.0)
        Vref = Blocks.Constant(k=1.0)
        Qref = Blocks.Constant(k=0.0)
        Qbranch = Blocks.Constant(k=0.0)
        Pbranch = Blocks.Constant(k=1.0)
        Vreg = Blocks.Constant(k=1.0)
        Ibranch = Blocks.Constant(k=1.0)
        Pe = Blocks.Constant(k=1.0)
        Pfaref = Blocks.Constant(k=1.0)
        Qgen = Blocks.Constant(k=0.0)
        Qgen0 = Blocks.Constant(k=0.0)
        Vdiff = Blocks.Constant(k=0)
    end
    @variables begin
        V_t(t), [description="Raw terminal voltage"]
        δ_v(t), [guess=0, description="voltage angle"]
        pir(t), [description="negative terminal current real part"]
        pii(t), [description="negative terminal current im part"]
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
        #Current injection from converter to terminal (negative for generation)
        [converter_interface.Ipout.u, converter_interface.Iqout.u].~ -1/CoB*[cos(δ_v) sin(δ_v); -sin(δ_v) cos(δ_v)] * [pir, pii]
        ##terminal.i_r ~ converter_interface.Ipout.u
        ##terminal.i_i ~ converter_interface.Iqout.u
        connect(plant_control.V_diff, Vdiff.output)
        connect(plant_control.freq, f.output)
        connect(plant_control.V_ref, Vref.output)
        connect(plant_control.Q_ref, Qref.output)
        connect(plant_control.Q_branch, Qbranch.output)
        connect(plant_control.P_branch, Pbranch.output)
        connect(plant_control.V_reg, Vreg.output)
        connect(plant_control.I_branch, Ibranch.output)
        connect(electrical_control.P_e, Pe.output)
        connect(electrical_control.P_faref, Pfaref.output)
        connect(electrical_control.Q_gen, Qgen.output)
        connect(converter_interface.Q_gen0, Qgen0.output)
        connect(plant_control.Pref_out, electrical_control.Pref_in)
        connect(plant_control.Qext_out, electrical_control.Qext_in)
        #electrical_control.Pref_in ~ 1
        #electrical_control.Qext_in ~ 1
        connect(electrical_control.Iqcmd_out, converter_interface.Iqcmd_in)
        connect(electrical_control.Ipcmd_out, converter_interface.Ipcmd_in)
        #converter_interface.Iqcmd_in.u ~ 1
        #converter_interface.Ipcmd_in.u ~ 1
    end
end