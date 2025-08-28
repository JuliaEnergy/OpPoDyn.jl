# Corrected WECC Large PV Plant Model for OpPoDyn.jl
# This fixes the SymbolicUtils error by properly structuring the model

using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using OpPoDyn

# First, define simplified component models that properly implement the interfaces
# You'll need to replace these with your actual Library.regc_a, Library.reec_b, Library.repc_a models

@mtkmodel regc_a begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        I_qrmax = 999.9
        I_qrmin = -999.9
        T_g = 0.02
        T_fltr = 0.01
        Brkpt = 0.9
        Zerox = 0.4
        L_vpl1 = 1.1
        L_vplsw = 0
        rrpwr = 10
    end
    @variables begin
        Vt_in(t), [input=true, description="Terminal voltage input"]
        Iqcmd_in(t), [input=true, description="Reactive current command input"]
        Ipcmd_in(t), [input=true, description="Active current command input"]
        Q_gen0(t), [input=true, description="Generator reactive power"]
    end
    @equations begin
        # Simplified converter interface equations
        # Replace with actual regc_a dynamics
        terminal.i_r ~ Ipcmd_in * cos(atan(terminal.u_i, terminal.u_r))
        terminal.i_i ~ Iqcmd_in * sin(atan(terminal.u_i, terminal.u_r))
    end
end

@mtkmodel reec_b begin
    @parameters begin
        V_dip = 0.0
        V_up = 1.1
        T_rv = 0.01
        V_ref0 = 0.95
        dbd1 = -0.01
        dbd2 = 0.01
        K_qv = 0.0
        I_qhl = 1.1
        I_qll = -1.1
        T_p = 0.01
        PfFlag = 0
        Q_min = -0.5
        Q_max = 0.5
        V_min = 0.85
        V_max = 1.15
        K_qp = 0.1
        K_qi = 0.02
        Vflag = 0
        K_vp = 0.1
        K_vi = 0.02
        QFlag = 0
        I_max = 1.0
        PqFlag = 0
        T_iq = 0.01
        T_pord = 0.02
        P_min = 0.0
        P_max = 1.0
        dP_min = -0.1
        dP_max = 0.1
    end
    @variables begin
        Vt_in(t), [input=true, description="Terminal voltage input"]
        Pref_in(t), [input=true, description="Active power reference input"]
        Qext_in(t), [input=true, description="External reactive power input"]
        Iqcmd_out(t), [output=true, description="Reactive current command output"]
        Ipcmd_out(t), [output=true, description="Active current command output"]
        P_e(t), [input=true, description="Electrical power"]
        P_faref(t), [input=true, description="Frequency ref"]
        Q_gen(t), [input=true, description="Generator reactive power"]
    end
    @equations begin
        # Simplified electrical control equations
        # Replace with actual reec_b dynamics
        Iqcmd_out ~ K_qp * (V_ref0 - Vt_in) + Qext_in
        Ipcmd_out ~ Pref_in
    end
end

@mtkmodel repc_a begin
    @parameters begin
        K_p = 0.1
        K_i = 0.02
        T_fltr = 0.01
        T_ft = 0.02
        T_fv = 0.02
        RefFlag = 0
        V_frz = 0.7
        R_c = 0.0
        X_c = 0.05
        K_c = 0.0
        VcombFlag = 0
        e_max = 0.1
        e_min = -0.1
        dbd_up = 0.0
        dbd_dn = 0.0
        Q_max = 0.5
        Q_min = -0.5
        K_pg = 0.1
        K_ig = 0.02
        T_p = 0.02
        fdbd1 = -0.02
        fdbd2 = 0.02
        femax = 0.1
        femin = -0.1
        P_max = 1.0
        P_min = 0.0
        T_lag = 0.02
        D_dn = 20.0
        D_up = 0.0
        P_plantref = 1.0
        freq_ref = 1.0
        freqFlag = 0
    end
    @components begin
        freq = Blocks.RealInput()
        V_ref = Blocks.RealInput()
        Q_ref = Blocks.RealInput()
        Q_branch = Blocks.RealInput()
        P_branch = Blocks.RealInput()
        V_reg = Blocks.RealInput()
        I_branch = Blocks.RealInput()
        Pref_out = Blocks.RealOutput()
        Qext_out = Blocks.RealOutput()
    end
    @equations begin
        # Simplified plant control equations
        # Replace with actual repc_a dynamics
        Pref_out.u ~ P_plantref
        Qext_out.u ~ K_p * (V_ref.u - 1.0)  # Simplified voltage control
    end
end

# Corrected WECC Large PV model
@mtkmodel WECC_large_PV begin
    @components begin
        terminal = Terminal()
        converter_interface = regc_a(
            I_qrmax = 999.9,
            I_qrmin = -999.9,
            T_g = 0.02,
            T_fltr = 0.01,
            Brkpt = 0.9,
            Zerox = 0.4,
            L_vpl1 = 1.1,
            L_vplsw = 0,
            rrpwr = 10)
        electrical_control = reec_b(
            V_dip = 0.0,
            V_up = 1.1,
            T_rv = 0.01,
            V_ref0 = 0.95,
            dbd1 = -0.01,
            dbd2 = 0.01,
            K_qv = 0.0,
            I_qhl = 1.1,
            I_qll = -1.1,
            T_p = 0.01,
            PfFlag = 0,
            Q_min = -0.5,
            Q_max = 0.5,
            V_min = 0.85,
            V_max = 1.15,
            K_qp = 0.1,
            K_qi = 0.02,
            Vflag = 0,
            K_vp = 0.1,
            K_vi = 0.02,
            QFlag = 0,
            I_max = 1.0,
            PqFlag = 0,
            T_iq = 0.01,
            T_pord = 0.02,
            P_min = 0.0,
            P_max = 1.0,
            dP_min = -0.1,
            dP_max = 0.1)
        plant_control = repc_a(
            K_p = 0.1,
            K_i = 0.02,
            T_fltr = 0.01,
            T_ft = 0.02,
            T_fv = 0.02,
            RefFlag = 0,
            V_frz = 0.7,
            R_c = 0.0,
            X_c = 0.05,
            K_c = 0.0,
            VcombFlag = 0,
            e_max = 0.1,
            e_min = -0.1,
            dbd_up = 0.0,
            dbd_dn = 0.0,
            Q_max = 0.5,
            Q_min = -0.5,
            K_pg = 0.1,
            K_ig = 0.02,
            T_p = 0.02,
            fdbd1 = -0.02,
            fdbd2 = 0.02,
            femax = 0.1,
            femin = -0.1,
            P_max = 1.0,
            P_min = 0.0,
            T_lag = 0.02,
            D_dn = 20.0,
            D_up = 0.0,
            P_plantref = 1.0,
            freq_ref = 1.0,
            freqFlag = 0)
        # Input constants
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
        # Calculate terminal voltage magnitude
        V_t ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        
        # Connect voltage to components (using equation, not connect for scalar signals)
        converter_interface.Vt_in ~ V_t
        electrical_control.Vt_in ~ V_t
        
        # Connect all the input signals using connect() for proper signal flow
        connect(plant_control.freq, f.output)
        connect(plant_control.V_ref, Vref.output)
        connect(plant_control.Q_ref, Qref.output)
        connect(plant_control.Q_branch, Qbranch.output)
        connect(plant_control.P_branch, Pbranch.output)
        connect(plant_control.V_reg, Vreg.output)
        connect(plant_control.I_branch, Ibranch.output)
        
        # Connect scalar inputs using equation (not connect for scalar variables)
        electrical_control.P_e ~ Pe.output.u
        electrical_control.P_faref ~ Pfaref.output.u
        electrical_control.Q_gen ~ Qgen.output.u
        converter_interface.Q_gen0 ~ Qgen0.output.u
        
        # Connect control chain outputs to inputs
        connect(plant_control.Pref_out, electrical_control.Pref_in)
        connect(plant_control.Qext_out, electrical_control.Qext_in)
        
        # Connect electrical control outputs to converter inputs
        electrical_control.Iqcmd_out ~ converter_interface.Iqcmd_in
        electrical_control.Ipcmd_out ~ converter_interface.Ipcmd_in
        
        # CRITICAL: Connect the main terminal to the converter terminal
        connect(terminal, converter_interface.terminal)
    end
end

# Corrected PVBus model
@mtkmodel PVBus begin
    @components begin
        busbar = BusBar()
        machine = WECC_large_PV()
    end
    @equations begin
        connect(busbar.terminal, machine.terminal)
    end
end

# Usage example:
# @named mtkbus2 = PVBus()
# bus_vertex = Bus(mtkbus2)