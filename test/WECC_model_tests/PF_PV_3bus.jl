EXPORT_FIGURES = true

using PowerDynamics
using OpPoDyn
using OpPoDyn.Library

using PowerDynamics.Library
using ModelingToolkit
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve

using CSV
using DataFrames
using CairoMakie
using Test


ref_pv = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV_pf","variables-testcase3Bus-with-event_addVariableNames.csv"),
    DataFrame;
    header = 3,
    decimal = ',',
    drop = (i, name) -> contains(string(name), "nrows="),
    silencewarnings = true
)


# bus 1 is provided from outside
PV_BUS = let
    ω_b = 2π*50

    # Powerflow results
    v_0 = 1.001047
    P_0 = 0.8888
    Q_0 = -0.3333

    @named PV = OpPoDyn.Library.WECC_large_PV_pf()
    busmodel = compile_bus(MTKBus(PV); current_source=true)
    compile_bus(busmodel, pf=pfPQ(P=P_0, Q=Q_0; current_source=true))
end

sol_pv = OpenIPSL_RePSSE_pv_pf_3bus(PV_BUS; ω_b = 2π*50);


## perform tests for all variables of interest
# Plant controls (repc_a)
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊repca₊P_ref), "repc_Pref") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊repca₊Q_ext), "repc_Qext") < 1e-3

# Electrical control (reec_b)
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊Q_measure), "repc_Qbranch") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊P_measure), "repc_Pbranch") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊V_t), "reec_Vt") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_pcmd), "reec_Ipcmd") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_qcmd), "reec_Iqcmd") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_pmax), "reec_Ipmax0") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_qmax), "reec_Iqmax0") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_qmin), "reec_Iqmin0") < 1e-3

# Renewable generator (regc_a)
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊regca₊I_lvpl), "regc_Ilvpl") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊I_measure), "repc_Ibranch") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊pvi), "repc_vregi") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊pvr), "repc_vregr") < 1e-3


# Create comprehensive comparison plot
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig = let
        fig = Figure(resolution=(1400, 1200))
        #ts = refine_timeseries(sol_pv.t)
        ts_refined = refine_timeseries(sol_pv.t, 20)
        ts = ts_refined[(ts_refined .>= 0) .& (ts_refined .<= 0.3)]
        idx = (ref_pv.time .>= 0) .& (ref_pv.time .<= 0.3)
        ts_filtered = ref_pv.time[idx]
        y_filtered  = ref_pv[!, Symbol("repc_Ibranch")][idx]

        # Plot 1: Current
        ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="[pu]", title="PV Generator Current State")
        lines!(ax1, ts_filtered,y_filtered; label="PowerFactory pir", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊I_measure)).u; label="PowerDynamics pir", color=:blue, linestyle=:dash, linewidth=2)
        axislegend(ax1)

        # Plot 2: pvi & pvr
        ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="[pu]", title="Generator States: pvi & pvr")
        lines!(ax2, ref_pv.time, ref_pv[!, Symbol("repc_vregi")]; label="PowerFactory pvi", color=:green, linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pvi)).u; label="PowerDynamics pvi", color=:green, linestyle=:dash, linewidth=2)
        lines!(ax2, ref_pv.time, ref_pv[!, Symbol("repc_vregr")]; label="PowerFactory pvr", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pvr)).u; label="PowerDynamics pvr", color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax2)

        # Plot 3: Vt
        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="Vt [pu]", title="Terminal Voltage Vt_in")
        lines!(ax3, ref_pv.time, ref_pv[!, Symbol("reec_Vt")]; label="PowerFactory Vt_in", color=:purple, linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊V_t)).u; label="PowerDynamics Vt_in", color=:purple, linestyle=:dash, linewidth=2)
        axislegend(ax3)

        # Plot 4: P_gen
        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="P [pu]", title="Generated Power P_gen")
        lines!(ax4, ref_pv.time, ref_pv[!, Symbol("repc_Pbranch")]; label="PowerFactory P_gen", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊P_measure)).u; label="PowerDynamics P_gen", color=:blue, linestyle=:dash, linewidth=2)
        axislegend(ax4)

        # Plot 5: Q_gen
        ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="Q [pu]", title="Generated Reactive Power Q_gen")
        lines!(ax5, ref_pv.time, ref_pv[!, Symbol("repc_Qbranch")]; label="PowerFactory Q_gen", color=:red, linewidth=2, alpha=0.7)
        lines!(ax5, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊Q_measure)).u; label="PowerDynamics Q_gen", color=:red, linestyle=:dash, linewidth=2)
        axislegend(ax5)

        # Plot 6: Ipcmd
        ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="Current [pu]", title="Ipcmd")
        lines!(ax6, ref_pv.time, ref_pv[!, Symbol("reec_Ipcmd")]; label="PowerFactory Ipcmd", color=:green, linewidth=2, alpha=0.7)
        lines!(ax6, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊reecb₊I_pcmd)).u; label="PowerDynamics Ipcmd", color=:green, linestyle=:dash, linewidth=2)
        axislegend(ax6)

        # Plot 7: Iqcmd
        ax7 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="Current [pu]", title="Iqcmd")
        lines!(ax7, ref_pv.time, ref_pv[!, Symbol("reec_Iqcmd")]; label="PowerFactory Iqcmd", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax7, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊reecb₊I_qcmd)).u; label="PowerDynamics Iqcmd", color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax7)

        # Plot 8: Qext & Pref (PlantController)
        ax8 = Axis(fig[4,2]; xlabel="Time [s]", ylabel="[pu]", title="PlantController: Qext & Pref")
        lines!(ax8, ref_pv.time, ref_pv[!, Symbol("repc_Qext")]; label="PowerFactory Qext", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax8, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊Q_ext)).u; label="PowerDynamics Qext", color=:blue, linestyle=:dash, linewidth=2)
        lines!(ax8, ref_pv.time, ref_pv[!, Symbol("repc_Pref")]; label="PowerFactory Pref", color=:red, linewidth=2, alpha=0.7)
        lines!(ax8, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊P_ref)).u; label="PowerDynamics Pref", color=:red, linestyle=:dash, linewidth=2)
        axislegend(ax8)

        fig
    end
    save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","PowerFactory_valid","PV_3bus_comparison_withLinefault_error.png"), fig)
end

# check Voltage dip and limits
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig_vdip = let
        ts_refined = refine_timeseries(sol_pv.t, 20)
        ts = ts_refined[(ts_refined .>= 0) .& (ts_refined .<= 0.3)]

        fig = Figure(resolution=(700, 400))
        ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="[-]", title="Voltage Dip Flag")
        lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊reecb₊Voltage_dip)).u; label="reec_b Voltage_dip", color=:blue, linewidth=2)
        lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊Voltage_dip)).u; label="repc_a Voltage_dip", color=:red, linewidth=2)
        lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊Q_lim)).u; label="repc_a Q_lim", color=:green, linewidth=2)
        lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊reecb₊P_lim)).u; label="reec_b P_lim", color=:pink, linewidth=2)
        lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊regca₊I_qr)).u; label="regc_a I_qr", color=:turquoise, linewidth=2)
        lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊regca₊I_pr)).u; label="regc_a I_pr", color=:orange, linewidth=2)
        #axislegend(ax)
        fig
    end
    save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","PowerFactory_valid","PV_3bus_voltage_dip.png"), fig_vdip)
end

# check different variables
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig_repc_filters = let
        ts_refined = refine_timeseries(sol_pv.t, 20)
        ts = ts_refined[(ts_refined .>= 0) .& (ts_refined .<= 0.3)]

        ymin = -1.75   # adjust as needed
        ymax =  1.2   # adjust as needed

        fig = Figure(resolution=(700, 400))
        ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="[pu]", title="repc_a: filter states vs. references")
        #lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊V_fltr)).u; label="V_fltr",  color=:pink,  linewidth=2)
        #lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊V_ref)).u;         label="V_ref",   color=:green,  linewidth=2)
        #lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊Q_fltr)).u;  label="Q_fltr",  color=:red,   linewidth=2)
        #lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊Q_ref)).u;          label="Q_ref",   color=:blue,   linewidth=2)
        lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊Q_measure)).u;  label="Q_measure / Q_branch",   color=:turquoise,   linewidth=2)
        ylims!(ax, ymin, ymax)
        axislegend(ax)
        fig
    end
    save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","PowerFactory_valid","PV_3bus_repc_filters.png"), fig_repc_filters)
end


# Input measures comparison: PowerFactory vs. simulation
ref_inputs = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV_pf","testcase3Bus-with-event_InputMeasures.csv"),
    DataFrame;
    header = 3,
    decimal = ',',
    silencewarnings = true
)
rename!(ref_inputs, 1 => :time)

if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig_inputs = let
        ts_refined = refine_timeseries(sol_pv.t, 20)
        ts = ts_refined[(ts_refined .>= 0) .& (ts_refined .<= 0.3)]
        idx = (ref_inputs.time .>= 0) .& (ref_inputs.time .<= 0.3)
        t_ref = ref_inputs.time[idx]

        fig = Figure(resolution=(1000, 800))

        # Plot 1: Frequency
        ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="[pu]", title="Frequency")
        lines!(ax1, t_ref, ref_inputs.f[idx];                                              label="PowerFactory f",   color=:blue, linewidth=2, alpha=0.7)
        lines!(ax1, ts,    sol_pv(ts, idxs=VIndex(:GEN1, :PV₊f)).u;                       label="PowerDynamics f",  color=:blue, linestyle=:dash, linewidth=2)
        axislegend(ax1)

        # Plot 2: Active and reactive power
        ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="[pu]", title="Active and Reactive Power")
        lines!(ax2, t_ref, ref_inputs.P_measure[idx];                                      label="PowerFactory P",   color=:blue,  linewidth=2, alpha=0.7)
        lines!(ax2, ts,    sol_pv(ts, idxs=VIndex(:GEN1, :PV₊P_measure)).u;               label="PowerDynamics P",  color=:blue,  linestyle=:dash, linewidth=2)
        lines!(ax2, t_ref, ref_inputs.Q_measure[idx];                                      label="PowerFactory Q",   color=:red,   linewidth=2, alpha=0.7)
        lines!(ax2, ts,    sol_pv(ts, idxs=VIndex(:GEN1, :PV₊Q_measure)).u;               label="PowerDynamics Q",  color=:red,   linestyle=:dash, linewidth=2)
        axislegend(ax2)

        # Plot 3: Currents
        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="[pu]", title="Currents")
        lines!(ax3, t_ref, ref_inputs.I_measure[idx];                                      label="PowerFactory |I|", color=:blue,   linewidth=2, alpha=0.7)
        lines!(ax3, ts,    sol_pv(ts, idxs=VIndex(:GEN1, :PV₊I_measure)).u;               label="PowerDynamics |I|",color=:blue,   linestyle=:dash, linewidth=2)
        lines!(ax3, t_ref, ref_inputs.pir[idx];                                            label="PowerFactory pir", color=:green,  linewidth=2, alpha=0.7)
        lines!(ax3, ts,    sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pir)).u;                     label="PowerDynamics pir",color=:green,  linestyle=:dash, linewidth=2)
        lines!(ax3, t_ref, ref_inputs.pii[idx];                                            label="PowerFactory pii", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax3, ts,    sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pii)).u;                     label="PowerDynamics pii",color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax3)

        # Plot 4: Voltages
        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="[pu]", title="Voltages")
        lines!(ax4, t_ref, ref_inputs.V_t[idx];                                            label="PowerFactory |V|", color=:blue,   linewidth=2, alpha=0.7)
        lines!(ax4, ts,    sol_pv(ts, idxs=VIndex(:GEN1, :PV₊V_t)).u;                     label="PowerDynamics |V|",color=:blue,   linestyle=:dash, linewidth=2)
        lines!(ax4, t_ref, ref_inputs.pvr[idx];                                            label="PowerFactory pvr", color=:green,  linewidth=2, alpha=0.7)
        lines!(ax4, ts,    sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pvr)).u;                     label="PowerDynamics pvr",color=:green,  linestyle=:dash, linewidth=2)
        lines!(ax4, t_ref, ref_inputs.pvi[idx];                                            label="PowerFactory pvi", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax4, ts,    sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pvi)).u;                     label="PowerDynamics pvi",color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax4)

        fig
    end
    save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","PowerFactory_valid","PV_3bus_input_measures.png"), fig_inputs)
end