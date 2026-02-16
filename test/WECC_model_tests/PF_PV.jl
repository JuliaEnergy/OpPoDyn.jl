using PowerDynamics
#PowerDynamics.load_pdtesting()
#using Main.PowerDynamicsTesting
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
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV","modelica_results_extended.csv"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# bus 1 is provided from outside
PV_BUS = let
    ω_b = 2π*50

    # Powerflow results
    v_0 = 1.0
    angle_0 = 0.0004339 #deg2rad(1.4753617387995086)
    P_0 = 0.015
    Q_0 = -0.056658

    @named PV = OpPoDyn.Library.WECC_large_PV_pf()
    busmodel = MTKBus(PV; name=:GEN1)
    #compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
    compile_bus(busmodel, pf=pfPV(V=v_0, P=P_0))
    #compile_bus(busmodel, pf=pfPQ(P=P_0, Q=Q_0))
end

sol_pv = OpenIPSL_RePSSE_pv(PV_BUS; ω_b = 2π*50);


## perform tests for all variables of interest
# Plant controls (repc_a)
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊repca₊P_ref), "pV.PlantController.Pref") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊repca₊Q_ext), "pV.PlantController.Qext") < 1e-3

# Electrical control (reec_b)
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊Q_gen), "pV.RenewableController.Qgen") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊P_gen), "pV.RenewableController.Pe") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊V_t), "pV.RenewableController.Vt") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_pcmd), "pV.RenewableController.Ipcmd") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_qcmd), "pV.RenewableController.Iqcmd") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_pmax), "pV.RenewableController.IPMAX.y") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_pmin), "pV.RenewableController.IPMIN.y") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_qmax), "pV.RenewableController.IQMAX.y") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊reecb₊I_qmin), "pV.RenewableController.IQMIN.y") < 1e-3

# Renewable generator (regc_a)
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊regca₊I_lvpl), "pV.RenewableGenerator.LVPL.y") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊pii), "pV.RenewableGenerator.p.ii") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊pir), "pV.RenewableGenerator.p.ir") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊pvi), "pV.RenewableGenerator.p.vi") < 1e-3
@test ref_rms_error(sol_pv, ref_pv, VIndex(:GEN1, :PV₊pvr), "pV.RenewableGenerator.p.vr") < 1e-3


# Create comprehensive comparison plot
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig = let
        fig = Figure(resolution=(1400, 1200))
        ts = refine_timeseries(sol_pv.t)

        # Plot 1: pir & pii
        ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="[pu]", title="PV Generator States: pir & pii")
        lines!(ax1, ref_pv.time, ref_pv[!, Symbol("pV.RenewableGenerator.p.ir")]; label="OpenIPSL pir", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pir)).u; label="PowerDynamics pir", color=:blue, linestyle=:dash, linewidth=2)
        lines!(ax1, ref_pv.time, ref_pv[!, Symbol("pV.RenewableGenerator.p.ii")]; label="OpenIPSL pii", color=:red, linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pii)).u; label="PowerDynamics pii", color=:red, linestyle=:dash, linewidth=2)
        axislegend(ax1)

        # Plot 2: pvi & pvr
        ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="[pu]", title="Generator States: pvi & pvr")
        lines!(ax2, ref_pv.time, ref_pv[!, Symbol("pV.RenewableGenerator.p.vi")]; label="OpenIPSL pvi", color=:green, linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pvi)).u; label="PowerDynamics pvi", color=:green, linestyle=:dash, linewidth=2)
        lines!(ax2, ref_pv.time, ref_pv[!, Symbol("pV.RenewableGenerator.p.vr")]; label="OpenIPSL pvr", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pvr)).u; label="PowerDynamics pvr", color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax2)

        # Plot 3: Vt
        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="Vt [pu]", title="Terminal Voltage Vt_in")
        lines!(ax3, ref_pv.time, ref_pv[!, Symbol("pV.RenewableController.Vt")]; label="OpenIPSL Vt_in", color=:purple, linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊V_t)).u; label="PowerDynamics Vt_in", color=:purple, linestyle=:dash, linewidth=2)
        axislegend(ax3)

        # Plot 4: P_gen
        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="P [pu]", title="Generated Power P_gen")
        lines!(ax4, ref_pv.time, ref_pv[!, Symbol("pV.RenewableController.Pe")]; label="OpenIPSL P_gen", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊P_gen)).u; label="PowerDynamics P_gen", color=:blue, linestyle=:dash, linewidth=2)
        axislegend(ax4)

        # Plot 5: Q_gen
        ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="Q [pu]", title="Generated Reactive Power Q_gen")
        lines!(ax5, ref_pv.time, ref_pv[!, Symbol("pV.RenewableController.Qgen")]; label="OpenIPSL Q_gen", color=:red, linewidth=2, alpha=0.7)
        lines!(ax5, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊Q_gen)).u; label="PowerDynamics Q_gen", color=:red, linestyle=:dash, linewidth=2)
        axislegend(ax5)

        # Plot 6: Ipcmd
        ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="Current [pu]", title="Ipcmd")
        lines!(ax6, ref_pv.time, ref_pv[!, Symbol("pV.RenewableController.Ipcmd")]; label="OpenIPSL Ipcmd", color=:green, linewidth=2, alpha=0.7)
        lines!(ax6, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊reecb₊I_pcmd)).u; label="PowerDynamics Ipcmd", color=:green, linestyle=:dash, linewidth=2)
        axislegend(ax6)

        # Plot 7: Iqcmd
        ax7 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="Current [pu]", title="Iqcmd")
        lines!(ax7, ref_pv.time, ref_pv[!, Symbol("pV.RenewableController.Iqcmd")]; label="OpenIPSL Iqcmd", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax7, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊reecb₊I_qcmd)).u; label="PowerDynamics Iqcmd", color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax7)

        # Plot 8: Qext & Pref (PlantController)
        ax8 = Axis(fig[4,2]; xlabel="Time [s]", ylabel="[pu]", title="PlantController: Qext & Pref")
        lines!(ax8, ref_pv.time, ref_pv[!, Symbol("pV.PlantController.Qext")]; label="OpenIPSL Qext", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax8, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊Q_ext)).u; label="PowerDynamics Qext", color=:blue, linestyle=:dash, linewidth=2)
        lines!(ax8, ref_pv.time, ref_pv[!, Symbol("pV.PlantController.Pref")]; label="OpenIPSL Pref", color=:red, linewidth=2, alpha=0.7)
        lines!(ax8, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊P_ref)).u; label="PowerDynamics Pref", color=:red, linestyle=:dash, linewidth=2)
        axislegend(ax8)

        fig
    end
    save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","OpenIPSL_valid","PV_comparison.png"), fig)
end


# --- Load extended reference data ---
ref_pv_extended = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV","modelica_results_extended.csv"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# --- Refine simulation timeseries ---
ts = refine_timeseries(sol_pv.t)

# --- PIR & PII ---
fig_pi = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="pu", title="PIR & PII Comparison")
    lines!(ax, ref_pv_extended.time, ref_pv_extended[!, "pV.RenewableGenerator.p.ir"]; label="OpenIPSL PIR", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pir)).u; label="PD PIR", color=Cycled(1), linewidth=2, linestyle=:dash)
    lines!(ax, ref_pv_extended.time, ref_pv_extended[!, "pV.RenewableGenerator.p.ii"]; label="OpenIPSL PII", color=Cycled(2), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pii)).u; label="PD PII", color=Cycled(2), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- PVI & PVR ---
fig_pv = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="pu", title="PVI & PVR Comparison")
    lines!(ax, ref_pv_extended.time, ref_pv_extended[!, "pV.RenewableGenerator.p.vi"]; label="OpenIPSL PVI", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pvi)).u; label="PD PVI", color=Cycled(1), linewidth=2, linestyle=:dash)
    lines!(ax, ref_pv_extended.time, ref_pv_extended[!, "pV.RenewableGenerator.p.vr"]; label="OpenIPSL PVR", color=Cycled(2), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊pvr)).u; label="PD PVR", color=Cycled(2), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- Terminal Voltage Vt_in.u ---
fig_Vt = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="Vt [pu]", title="Terminal Voltage Vt Comparison")
    lines!(ax, ref_pv_extended.time, ref_pv_extended[!, "pV.RenewableController.Vt"]; label="OpenIPSL Vt", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊V_t)).u; label="PD Vt", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- Active Power P_gen ---
fig_Pgen = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="P [pu]", title="Active Power P_gen Comparison")
    lines!(ax, ref_pv_extended.time, ref_pv_extended[!, "pV.RenewableController.Pe"]; label="OpenIPSL P_gen", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊P_gen)).u; label="PD P_gen", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- Reactive Power Q_gen ---
fig_Qgen = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="Q [pu]", title="Reactive Power Q_gen Comparison")
    lines!(ax, ref_pv_extended.time, ref_pv_extended[!, "pV.RenewableController.Qgen"]; label="OpenIPSL Q_gen", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊Q_gen)).u; label="PD Q_gen", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- Ipcmd ---
fig_Ipcmd = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="I [pu]", title="Ipcmd Comparison")
    lines!(ax, ref_pv_extended.time, ref_pv_extended[!, "pV.RenewableController.Ipcmd"]; label="OpenIPSL Ipcmd", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊reecb₊I_pcmd)).u; label="PD Ipcmd", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- Iqcmd ---
fig_Iqcmd = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="I [pu]", title="Iqcmd Comparison")
    lines!(ax, ref_pv_extended.time, ref_pv_extended[!, "pV.RenewableController.Iqcmd"]; label="OpenIPSL Iqcmd", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊reecb₊I_qcmd)).u; label="PD Iqcmd", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- PlantController Qext & Pref ---
fig_plant = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="pu", title="PlantController: Qext & Pref")
    lines!(ax, ref_pv_extended.time, ref_pv_extended[!, "pV.PlantController.Qext"]; label="OpenIPSL Qext", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊Q_ext)).u; label="PD Qext", color=Cycled(1), linewidth=2, linestyle=:dash)
    lines!(ax, ref_pv_extended.time, ref_pv_extended[!, "pV.PlantController.Pref"]; label="OpenIPSL Pref", color=Cycled(2), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊P_ref)).u; label="PD Pref", color=Cycled(2), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

#voltage V_dip
fig_plant = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="pu", title="Voltage_dip")
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊repca₊Voltage_dip)).u; label="PD Voltage_dip plant control", color=Cycled(2), linewidth=2, linestyle=:dash)
    lines!(ax, ts, sol_pv(ts, idxs=VIndex(:GEN1, :PV₊reecb₊Voltage_dip)).u; label="PD Voltage_dip electrical control", color=Cycled(2), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end
