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

ref_bess = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","BESS","modelica_results_extended.csv"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# bus 1 is provided from outside
BESS_BUS = let
    ω_b = 2π*50

    # Powerflow results
    v_0 = 1.0
    angle_0 = deg2rad(1.4753617387995086) #in rad
    P_0 = 0.015
    Q_0 = -0.056658

    @named BESS = OpPoDyn.Library.WECC_BESS()
    busmodel = MTKBus(BESS; name=:GEN1)
    #compile_bus(busmodel, pf=pfSlack(V=v_0, δ=angle_0))
    compile_bus(busmodel, pf=pfPQ(P=P_0, Q=Q_0))
end

sol_bess = OpenIPSL_RePSSE(BESS_BUS);

## perform tests for all variables of interest
# Plant controls (repc_a)
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊plant_control₊P_ref), "bESS.PlantController.Pref") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊plant_control₊Q_ext), "bESS.PlantController.Qext") < 1e-3

# Electrical control (reec_b)
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊Q_gen), "bESS.RenewableController.Qgen") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊P_gen), "bESS.RenewableController.Pe") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊V_t), "bESS.RenewableController.Vt") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊electrical_control₊I_pcmd), "bESS.RenewableController.Ipcmd") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊electrical_control₊I_qcmd), "bESS.RenewableController.Iqcmd") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊electrical_control₊soc_lim), "bESS.RenewableController.sOC_logic.SOC") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊electrical_control₊I_pmax), "bESS.RenewableController.IPMAX.y") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊electrical_control₊I_pmin), "bESS.RenewableController.IPMIN.y") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊electrical_control₊I_qmax), "bESS.RenewableController.IQMAX.y") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊electrical_control₊I_qmin), "bESS.RenewableController.IQMIN.y") < 1e-3

# Renewable generator (regc_a)
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊converter_interface₊I_lvpl), "bESS.RenewableGenerator.LVPL.y") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊pii), "bESS.RenewableGenerator.p.ii") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊pir), "bESS.RenewableGenerator.p.ir") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊pvi), "bESS.RenewableGenerator.p.vi") < 1e-3
@test ref_rms_error(sol_bess, ref_bess, VIndex(:GEN1, :BESS₊pvr), "bESS.RenewableGenerator.p.vr") < 1e-3


# Create comprehensive comparison plot
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig = let
        fig = Figure(resolution=(1400, 1200))
        ts = refine_timeseries(sol_bess.t)

        # Plot 1: pir & pii
        ax1 = Axis(fig[1,1]; xlabel="Time [s]", ylabel="[pu]", title="Generator States: pir & pii")
        lines!(ax1, ref_bess.time, ref_bess[!, Symbol("bESS.RenewableGenerator.p.ir")]; label="OpenIPSL pir", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊pir)).u; label="PowerDynamics pir", color=:blue, linestyle=:dash, linewidth=2)
        lines!(ax1, ref_bess.time, ref_bess[!, Symbol("bESS.RenewableGenerator.p.ii")]; label="OpenIPSL pii", color=:red, linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊pii)).u; label="PowerDynamics pii", color=:red, linestyle=:dash, linewidth=2)
        axislegend(ax1)

        # Plot 2: pvi & pvr
        ax2 = Axis(fig[1,2]; xlabel="Time [s]", ylabel="[pu]", title="Generator States: pvi & pvr")
        lines!(ax2, ref_bess.time, ref_bess[!, Symbol("bESS.RenewableGenerator.p.vi")]; label="OpenIPSL pvi", color=:green, linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊pvi)).u; label="PowerDynamics pvi", color=:green, linestyle=:dash, linewidth=2)
        lines!(ax2, ref_bess.time, ref_bess[!, Symbol("bESS.RenewableGenerator.p.vr")]; label="OpenIPSL pvr", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊pvr)).u; label="PowerDynamics pvr", color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax2)

        # Plot 3: Vt_in
        ax3 = Axis(fig[2,1]; xlabel="Time [s]", ylabel="Vt [pu]", title="Terminal Voltage Vt_in")
        lines!(ax3, ref_bess.time, ref_bess[!, Symbol("bESS.RenewableController.Vt")]; label="OpenIPSL Vt_in", color=:purple, linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊V_t)).u; label="PowerDynamics Vt_in", color=:purple, linestyle=:dash, linewidth=2)
        axislegend(ax3)

        # Plot 4: P_gen
        ax4 = Axis(fig[2,2]; xlabel="Time [s]", ylabel="P [pu]", title="Generated Power P_gen")
        lines!(ax4, ref_bess.time, ref_bess[!, Symbol("bESS.RenewableController.Pe")]; label="OpenIPSL P_gen", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊P_gen)).u; label="PowerDynamics P_gen", color=:blue, linestyle=:dash, linewidth=2)
        axislegend(ax4)

        # Plot 5: Q_gen
        ax5 = Axis(fig[3,1]; xlabel="Time [s]", ylabel="Q [pu]", title="Generated Reactive Power Q_gen")
        lines!(ax5, ref_bess.time, ref_bess[!, Symbol("bESS.RenewableController.Qgen")]; label="OpenIPSL Q_gen", color=:red, linewidth=2, alpha=0.7)
        lines!(ax5, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊Q_gen)).u; label="PowerDynamics Q_gen", color=:red, linestyle=:dash, linewidth=2)
        axislegend(ax5)

        # Plot 6: Ipcmd
        ax6 = Axis(fig[3,2]; xlabel="Time [s]", ylabel="Current [pu]", title="Ipcmd")
        lines!(ax6, ref_bess.time, ref_bess[!, Symbol("bESS.RenewableController.Ipcmd")]; label="OpenIPSL Ipcmd", color=:green, linewidth=2, alpha=0.7)
        lines!(ax6, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊electrical_control₊I_pcmd)).u; label="PowerDynamics Ipcmd", color=:green, linestyle=:dash, linewidth=2)
        axislegend(ax6)

        # Plot 7: Iqcmd
        ax7 = Axis(fig[4,1]; xlabel="Time [s]", ylabel="Current [pu]", title="Iqcmd")
        lines!(ax7, ref_bess.time, ref_bess[!, Symbol("bESS.RenewableController.Iqcmd")]; label="OpenIPSL Iqcmd", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax7, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊electrical_control₊I_qcmd)).u; label="PowerDynamics Iqcmd", color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax7)

        # Plot 8: Qext & Pref (PlantController)
        ax8 = Axis(fig[4,2]; xlabel="Time [s]", ylabel="[pu]", title="PlantController: Qext & Pref")
        lines!(ax8, ref_bess.time, ref_bess[!, Symbol("bESS.PlantController.Qext")]; label="OpenIPSL Qext", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax8, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊plant_control₊Q_ext)).u; label="PowerDynamics Qext", color=:blue, linestyle=:dash, linewidth=2)
        lines!(ax8, ref_bess.time, ref_bess[!, Symbol("bESS.PlantController.Pref")]; label="OpenIPSL Pref", color=:red, linewidth=2, alpha=0.7)
        lines!(ax8, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊plant_control₊P_ref)).u; label="PowerDynamics Pref", color=:red, linestyle=:dash, linewidth=2)
        axislegend(ax8)

        fig
    end
    save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","OpenIPSL_valid","BESS_comparison.png"), fig)
end


# --- Load extended reference data ---
ref_bess_extended = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","BESS","modelica_results_extended.csv"),
    DataFrame;
    drop=(i,name) -> contains(string(name), "nrows="),
    silencewarnings=true
)

# --- Refine simulation timeseries ---
ts = refine_timeseries(sol_bess.t)

# --- PIR & PII ---
fig_pi = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="pu", title="PIR & PII Comparison")
    lines!(ax, ref_bess_extended.time, ref_bess_extended[!, "bESS.RenewableGenerator.p.ir"]; label="OpenIPSL PIR", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊pir)).u; label="PD PIR", color=Cycled(1), linewidth=2, linestyle=:dash)
    lines!(ax, ref_bess_extended.time, ref_bess_extended[!, "bESS.RenewableGenerator.p.ii"]; label="OpenIPSL PII", color=Cycled(2), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊pii)).u; label="PD PII", color=Cycled(2), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- PVI & PVR ---
fig_pv = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="pu", title="PVI & PVR Comparison")
    lines!(ax, ref_bess_extended.time, ref_bess_extended[!, "bESS.RenewableGenerator.p.vi"]; label="OpenIPSL PVI", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊pvi)).u; label="PD PVI", color=Cycled(1), linewidth=2, linestyle=:dash)
    lines!(ax, ref_bess_extended.time, ref_bess_extended[!, "bESS.RenewableGenerator.p.vr"]; label="OpenIPSL PVR", color=Cycled(2), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊pvr)).u; label="PD PVR", color=Cycled(2), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- Terminal Voltage Vt_in.u ---
fig_Vt = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="Vt [pu]", title="Terminal Voltage Vt Comparison")
    lines!(ax, ref_bess_extended.time, ref_bess_extended[!, "bESS.RenewableController.Vt"]; label="OpenIPSL Vt", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊V_t)).u; label="PD Vt", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- Active Power P_gen ---
fig_Pgen = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="P [pu]", title="Active Power P_gen Comparison")
    lines!(ax, ref_bess_extended.time, ref_bess_extended[!, "bESS.RenewableController.Pe"]; label="OpenIPSL P_gen", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊P_gen)).u; label="PD P_gen", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- Reactive Power Q_gen ---
fig_Qgen = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="Q [pu]", title="Reactive Power Q_gen Comparison")
    lines!(ax, ref_bess_extended.time, ref_bess_extended[!, "bESS.RenewableController.Qgen"]; label="OpenIPSL Q_gen", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊Q_gen)).u; label="PD Q_gen", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- Ipcmd ---
fig_Ipcmd = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="I [pu]", title="Ipcmd Comparison")
    lines!(ax, ref_bess_extended.time, ref_bess_extended[!, "bESS.RenewableController.Ipcmd"]; label="OpenIPSL Ipcmd", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊electrical_control₊I_pcmd)).u; label="PD Ipcmd", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- Iqcmd ---
fig_Iqcmd = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="I [pu]", title="Iqcmd Comparison")
    lines!(ax, ref_bess_extended.time, ref_bess_extended[!, "bESS.RenewableController.Iqcmd"]; label="OpenIPSL Iqcmd", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊electrical_control₊I_qcmd)).u; label="PD Iqcmd", color=Cycled(1), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end

# --- PlantController Qext & Pref ---
fig_plant = let
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1]; xlabel="Time [s]", ylabel="pu", title="PlantController: Qext & Pref")
    lines!(ax, ref_bess_extended.time, ref_bess_extended[!, "bESS.PlantController.Qext"]; label="OpenIPSL Qext", color=Cycled(1), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊plant_control₊Q_ext)).u; label="PD Qext", color=Cycled(1), linewidth=2, linestyle=:dash)
    lines!(ax, ref_bess_extended.time, ref_bess_extended[!, "bESS.PlantController.Pref"]; label="OpenIPSL Pref", color=Cycled(2), linewidth=2, alpha=0.5)
    lines!(ax, ts, sol_bess(ts, idxs=VIndex(:GEN1, :BESS₊plant_control₊P_ref)).u; label="PD Pref", color=Cycled(2), linewidth=2, linestyle=:dash)
    axislegend(ax; position=:rt)
    fig
end
