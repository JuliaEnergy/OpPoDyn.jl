# WECC Model for  large scale PV plant taken from https://www.esig.energy/wiki-main-page/generic-models-pv-plants/#WECC_Generic_Model_for_Large-scale_PV_Plants
# default parameter taken from typcial values, flags=0

@mtkmodel WECC_large_PV begin
    @components begin
        terminal=Terminal()
        converter_interface = regc_a(
            I_qrmax = 999.9,
            I_qrmin = -999.9,
            T_g     = 0.02,
            T_fltr  = 0.01,
            Brkpt   = 0.9,
            Zerox   = 0.4,
            L_vpl1  = 1.1,
            L_vplsw = 0,
            rrpwr   = 10,
            V_0lim = 1.2,
            K_hv = 0.7)
        #=electrical_control = Library.reec_b(
            V_dip   = 0.0,
            V_up    = 1.1,
            T_rv    = 0.01,
            V_ref0  = 0.95,
            dbd1    = -0.01,
            dbd2    = 0.01,
            K_qv    = 0.0,
            I_qhl   = 1.1,
            I_qll   = -1.1,
            T_p     = 0.01,
            PfFlag  = 0,
            Q_min   = -0.5, #no examplary values given
            Q_max   = 0.5, #no examplary values given
            V_min   = 0.85,
            V_max   = 1.15,
            K_qp    = 0.1, #no examplary values given
            K_qi    = 0.02, #no examplary values given
            Vflag   = 0,
            K_vp    = 0.1, #no examplary values given
            K_vi    = 0.02, #no examplary values given
            QFlag   = 0,
            I_max   = 1.0,
            PqFlag  = 0,
            T_iq    = 0.01,
            T_pord  = 0.02, #no examplary values given
            P_min   = 0.0,
            P_max   = 1.0,
            dP_min  = -0.1, #no examplary values given
            dP_max  = 0.1) #no examplary values given
        plant_control = Library.repc_a(
            K_p        = 0.1, #no examplary values given
            K_i        = 0.02, #no examplary values given
            T_fltr     = 0.01,
            T_ft       = 0.02, #no examplary values given
            T_fv       = 0.02, #no examplary values given
            RefFlag    = 0,
            V_frz      = 0.7,
            R_c        = 0.0,
            X_c        = 0.05,
            K_c        = 0.0, #no examplary values given
            VcombFlag  = 0,
            e_max      = 0.1, #no examplary values given
            e_min      = -0.1, #no examplary values given
            dbd_up     = 0.0,
            dbd_dn     = 0.0,
            Q_max      = 0.5, #no examplary values given
            Q_min      = -0.5, #no examplary values given
            K_pg       = 0.1, #no examplary values given
            K_ig       = 0.02, #no examplary values given
            T_p        = 0.02, #no examplary values given
            fdbd1      = -0.02, #no examplary values given
            fdbd2      = 0.02, #no examplary values given
            femax      = 0.1, #no examplary values given
            femin      = -0.1, #no examplary values given
            P_max      = 1.0, #no examplary values given
            P_min      = 0.0, #no examplary values given
            T_lag      = 0.02, #no examplary values given
            D_dn       = 20.0,
            D_up       = 0.0,
            P_plantref = 1.0, #from powerflow
            freq_ref   = 1.0,
            freqFlag   = 0)=#
        f = Blocks.Constant(k=1.0)
        Vref = Blocks.Constant(k=1.0)
        Qref = Blocks.Constant(k=1.0)
        Qbranch = Blocks.Constant(k=1.0)
        Pbranch = Blocks.Constant(k=1.0)
        Vreg = Blocks.Constant(k=1.0)
        Ibranch = Blocks.Constant(k=1.0)
        Pe = Blocks.Constant(k=1.0)
        Pfaref = Blocks.Constant(k=1.0)
        Qgen = Blocks.Constant(k=1.0)
        Qgen0 = Blocks.Constant(k=1.0)
    end
    @variables begin
        V_t(t), [description="Raw terminal voltage"]
    end
    @equations begin
        V_t ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        converter_interface.Vt_in.u ~ V_t
        #electrical_control.Vt_in.u ~ V_t
        ##connect(converter_interface.Vt_in, V_t)
        ##connect(electrical_control.Vt_in, V_t)
        #Current injection from converter to terminal (negative for generation)
        #terminal.i_r ~ converter_interface.Ipout.u
        #terminal.i_i ~ -converter_interface.Iqout.u
        #connect(plant_control.freq, f.output)
        #connect(plant_control.V_ref, Vref.output)
        #connect(plant_control.Q_ref, Qref.output)
        #connect(plant_control.Q_branch, Qbranch.output)
        #connect(plant_control.P_branch, Pbranch.output)
        #connect(plant_control.V_reg, Vreg.output)
        #connect(plant_control.I_branch, Ibranch.output)
        #connect(electrical_control.P_e, Pe.output)
        #connect(electrical_control.P_faref, Pfaref.output)
        #connect(electrical_control.Q_gen, Qgen.output)
        connect(converter_interface.Q_gen0, Qgen0.output)
        #connect(plant_control.Pref_out, electrical_control.Pref_in)
        #connect(plant_control.Qext_out, electrical_control.Qext_in)
        #electrical_control.Pref_in ~ 1
        #electrical_control.Qext_in ~ 1
        #connect(electrical_control.Iqcmd_out, converter_interface.Iqcmd_in)
        #connect(electrical_control.Ipcmd_out, converter_interface.Ipcmd_in)
        converter_interface.Iqcmd_in.u ~ 1
        converter_interface.Ipcmd_in.u ~ 1
    end
end