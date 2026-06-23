EXPORT_FIGURES = true

# ── Pakete ────────────────────────────────────────────────────────────────────
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
using OrderedCollections
using Test


# ══════════════════════════════════════════════════════════════════════════════
# ABSCHNITT 1: Vollständige Netzwerksimulation (Standardfall)
#   Julia-Modell ohne Switch-off-Logik, verglichen mit PF ohne Switch-off-Schwellwert
# ══════════════════════════════════════════════════════════════════════════════

# Referenz interne Variablen:
#   PF-Vollsimulation, Switch-off-Schwellwert = 0 pu (deaktiviert), mit Linienfehler
ref_fullsim = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV_pf",
        "variables-testcase3Bus-with-event_addVariableNames_withoutThresholdOperationVoltage.csv"),
    DataFrame;
    header = 3,
    decimal = ',',
    drop = (i, name) -> contains(string(name), "nrows="),
    silencewarnings = true
)

# Referenz Klemmengrößen (Spannung, Strom, P, Q):
#   PF-Vollsimulation, Switch-off-Schwellwert = 0 pu (deaktiviert), mit Linienfehler
ref_terminal_nothresh = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV_pf",
        "testcase3Bus-with-event_InputMeasures_withoutThresholdOperationVoltage.csv"),
    DataFrame;
    header = 3,
    delim = ';',
    decimal = ',',
    silencewarnings = true
)
rename!(ref_terminal_nothresh, 1 => :time)

# ── Simulation ────────────────────────────────────────────────────────────────

# PV-Bus aufbauen (PowerFlow-Ergebnisse aus PF)
PV_BUS = let
    v_0 = 1.001047
    P_0 = 0.8888
    Q_0 = -0.3333
    @named PV = OpPoDyn.Library.WECC_large_PV_pf()
    busmodel = compile_bus(MTKBus(PV); current_source=true)
    compile_bus(busmodel, pf=pfPQ(P=P_0, Q=Q_0; current_source=true))
end

# Lösung: vollständige Netzwerksimulation, keine fixen Inputs, kein Switch-off
sol_fullsim = OpenIPSL_RePSSE_pv_pf_3bus(PV_BUS; ω_b = 2π*50)

# ── Tests ─────────────────────────────────────────────────────────────────────

# Plant controls (repc_a)
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊repca₊P_ref), "repc_Pref")  < 1e-3
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊repca₊Q_ext), "repc_Qext")  < 1e-3

# Electrical control (reec_b)
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊Q_measure),      "repc_Qbranch") < 1e-3
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊P_measure),      "repc_Pbranch") < 1e-3
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊V_t),            "reec_Vt")      < 1e-3
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊reecb₊I_pcmd),   "reec_Ipcmd")   < 1e-3
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊reecb₊I_qcmd),   "reec_Iqcmd")   < 1e-3
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊reecb₊I_pmax),   "reec_Ipmax0")  < 1e-3
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊reecb₊I_qmax),   "reec_Iqmax0")  < 1e-3
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊reecb₊I_qmin),   "reec_Iqmin0")  < 1e-3

# Renewable generator (regc_a)
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊regca₊I_lvpl), "regc_Ilvpl")   < 1e-3
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊I_measure),    "repc_Ibranch") < 1e-3
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊pvi),          "repc_vregi")   < 1e-3
@test ref_rms_error(sol_fullsim, ref_fullsim, VIndex(:GEN1, :PV₊pvr),          "repc_vregr")   < 1e-3

# ── Plots ─────────────────────────────────────────────────────────────────────

# Überblick: I_measure, pvr/pvi, Vt, P, Q, Ipcmd, Iqcmd, Qext/Pref
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig_fullsim_overview = let
        ts  = let r = refine_timeseries(sol_fullsim.t, 20); r[(r .>= 0) .& (r .<= 0.3)] end
        idx = (ref_fullsim.time .>= 0) .& (ref_fullsim.time .<= 0.3)

        fig = Figure(resolution=(1400, 1200))

        ax1 = Axis(fig[1,1]; xlabel="Zeit [s]", ylabel="|I| [pu]", title="Stromstärke am Einspeisepunkt |I_branch| = √(pir²+pii²)")
        lines!(ax1, ref_fullsim.time[idx], ref_fullsim.repc_Ibranch[idx]; label="PF  repc_Ibranch", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊I_measure)).u; label="Julia  I_measure", color=:blue, linestyle=:dash, linewidth=2)
        axislegend(ax1)

        ax2 = Axis(fig[1,2]; xlabel="Zeit [s]", ylabel="u [pu]", title="Regler-Referenzspannung: Realteil pvr und Imaginärteil pvi")
        lines!(ax2, ref_fullsim.time[idx], ref_fullsim.repc_vregi[idx]; label="PF  pvi = Im{v_reg}", color=:green,  linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊pvi)).u; label="Julia  pvi", color=:green,  linestyle=:dash, linewidth=2)
        lines!(ax2, ref_fullsim.time[idx], ref_fullsim.repc_vregr[idx]; label="PF  pvr = Re{v_reg}", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊pvr)).u; label="Julia  pvr", color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax2)

        ax3 = Axis(fig[2,1]; xlabel="Zeit [s]", ylabel="|V_t| [pu]", title="Klemmenspannungsbetrag |V_t| (Eingang REEC-B)")
        lines!(ax3, ref_fullsim.time[idx], ref_fullsim.reec_Vt[idx]; label="PF  reec_Vt", color=:purple, linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊V_t)).u; label="Julia  V_t", color=:purple, linestyle=:dash, linewidth=2)
        axislegend(ax3)

        ax4 = Axis(fig[2,2]; xlabel="Zeit [s]", ylabel="P [pu]", title="Wirkleistung P am Einspeisepunkt (repc_Pbranch)")
        lines!(ax4, ref_fullsim.time[idx], ref_fullsim.repc_Pbranch[idx]; label="PF  repc_Pbranch", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊P_measure)).u; label="Julia  P_measure", color=:blue, linestyle=:dash, linewidth=2)
        axislegend(ax4)

        ax5 = Axis(fig[3,1]; xlabel="Zeit [s]", ylabel="Q [pu]", title="Blindleistung Q am Einspeisepunkt (repc_Qbranch)")
        lines!(ax5, ref_fullsim.time[idx], ref_fullsim.repc_Qbranch[idx]; label="PF  repc_Qbranch", color=:red, linewidth=2, alpha=0.7)
        lines!(ax5, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊Q_measure)).u; label="Julia  Q_measure", color=:red, linestyle=:dash, linewidth=2)
        axislegend(ax5)

        ax6 = Axis(fig[3,2]; xlabel="Zeit [s]", ylabel="I_pcmd [pu]", title="Wirkstrom-Sollwert I_pcmd (Ausgang REEC-B → REGC-A)")
        lines!(ax6, ref_fullsim.time[idx], ref_fullsim.reec_Ipcmd[idx]; label="PF  reec_Ipcmd", color=:green, linewidth=2, alpha=0.7)
        lines!(ax6, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊reecb₊I_pcmd)).u; label="Julia  reecb.I_pcmd", color=:green, linestyle=:dash, linewidth=2)
        axislegend(ax6)

        ax7 = Axis(fig[4,1]; xlabel="Zeit [s]", ylabel="I_qcmd [pu]", title="Blindstrom-Sollwert I_qcmd (Ausgang REEC-B → REGC-A)")
        lines!(ax7, ref_fullsim.time[idx], ref_fullsim.reec_Iqcmd[idx]; label="PF  reec_Iqcmd", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax7, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊reecb₊I_qcmd)).u; label="Julia  reecb.I_qcmd", color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax7)

        ax8 = Axis(fig[4,2]; xlabel="Zeit [s]", ylabel="[pu]", title="Anlagenregler REPC-A: Blindleistungssollwert Q_ext und Wirkleistungsreferenz P_ref")
        lines!(ax8, ref_fullsim.time[idx], ref_fullsim.repc_Qext[idx]; label="PF  repc_Qext (REPC-A Ausgang)",  color=:blue, linewidth=2, alpha=0.7)
        lines!(ax8, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊repca₊Q_ext)).u; label="Julia  repca.Q_ext", color=:blue, linestyle=:dash, linewidth=2)
        lines!(ax8, ref_fullsim.time[idx], ref_fullsim.repc_Pref[idx]; label="PF  repc_Pref (REPC-A Ausgang)",  color=:red,  linewidth=2, alpha=0.7)
        lines!(ax8, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊repca₊P_ref)).u; label="Julia  repca.P_ref", color=:red,  linestyle=:dash, linewidth=2)
        axislegend(ax8)

        fig
    end
    save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","PowerFactory_valid","PV_3bus_fullsim_overview.png"), fig_fullsim_overview)
end

# Klemmengrößen: pir, pii, pvr, pvi, |I|, |V|, P, Q
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig_fullsim_terminal = let
        ts  = let r = refine_timeseries(sol_fullsim.t, 20); r[(r .>= 0) .& (r .<= 0.3)] end
        idx = (ref_terminal_nothresh.time .>= 0) .& (ref_terminal_nothresh.time .<= 0.3)
        t_ref = ref_terminal_nothresh.time[idx]

        fig = Figure(resolution=(1000, 800))

        ax1 = Axis(fig[1,1]; xlabel="Zeit [s]", ylabel="f [pu]", title="Frequenz am Anschlussknoten (Eingang REPC-A)")
        lines!(ax1, t_ref, ref_terminal_nothresh.f[idx]; label="PF  f", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax1, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊f)).u; label="Julia  f", color=:blue, linestyle=:dash, linewidth=2)
        axislegend(ax1)

        ax2 = Axis(fig[1,2]; xlabel="Zeit [s]", ylabel="[pu]", title="Wirk- und Blindleistung P, Q am Einspeisepunkt")
        lines!(ax2, t_ref, ref_terminal_nothresh.P_measure[idx]; label="PF  P_measure", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊P_measure)).u; label="Julia  P_measure", color=:blue, linestyle=:dash, linewidth=2)
        lines!(ax2, t_ref, ref_terminal_nothresh.Q_measure[idx]; label="PF  Q_measure", color=:red,  linewidth=2, alpha=0.7)
        lines!(ax2, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊Q_measure)).u; label="Julia  Q_measure", color=:red,  linestyle=:dash, linewidth=2)
        axislegend(ax2)

        ax3 = Axis(fig[2,1]; xlabel="Zeit [s]", ylabel="[pu]", title="Klemmenströme: Betrag |I| sowie Realteil pir und Imaginärteil pii")
        lines!(ax3, t_ref, ref_terminal_nothresh.I_measure[idx]; label="PF  |I_branch| = √(pir²+pii²)", color=:blue,   linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊I_measure)).u; label="Julia  I_measure", color=:blue,   linestyle=:dash, linewidth=2)
        lines!(ax3, t_ref, ref_terminal_nothresh.pir[idx]; label="PF  pir = Re{i}", color=:green,  linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊pir)).u; label="Julia  pir", color=:green,  linestyle=:dash, linewidth=2)
        lines!(ax3, t_ref, ref_terminal_nothresh.pii[idx]; label="PF  pii = Im{i}", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax3, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊pii)).u; label="Julia  pii", color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax3)

        ax4 = Axis(fig[2,2]; xlabel="Zeit [s]", ylabel="[pu]", title="Klemmenspannungen: Betrag |V_t| sowie Realteil pvr und Imaginärteil pvi")
        lines!(ax4, t_ref, ref_terminal_nothresh.V_t[idx]; label="PF  |V_t| = √(pvr²+pvi²)", color=:blue,   linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊V_t)).u; label="Julia  V_t", color=:blue,   linestyle=:dash, linewidth=2)
        lines!(ax4, t_ref, ref_terminal_nothresh.pvr[idx]; label="PF  pvr = Re{v}", color=:green,  linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊pvr)).u; label="Julia  pvr", color=:green,  linestyle=:dash, linewidth=2)
        lines!(ax4, t_ref, ref_terminal_nothresh.pvi[idx]; label="PF  pvi = Im{v}", color=:orange, linewidth=2, alpha=0.7)
        lines!(ax4, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊pvi)).u; label="Julia  pvi", color=:orange, linestyle=:dash, linewidth=2)
        axislegend(ax4)

        fig
    end
    save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","PowerFactory_valid","PV_3bus_fullsim_terminal.png"), fig_fullsim_terminal)
end

# Bus-A-Spannung (Knoten 1): Betrag, Winkel, Real- und Imaginärteil
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig_fullsim_busA = let
        ts  = let r = refine_timeseries(sol_fullsim.t, 20); r[(r .>= 0.05) .& (r .<= 0.5)] end
        idx   = (ref_fullsim.time .>= 0.05) .& (ref_fullsim.time .<= 0.5)
        t_ref = ref_fullsim.time[idx]
        ur_jl = sol_fullsim(ts, idxs=VIndex(1, :busbar₊u_r)).u
        ui_jl = sol_fullsim(ts, idxs=VIndex(1, :busbar₊u_i)).u

        fig = Figure(resolution=(1100, 900))

        ax1 = Axis(fig[1,1]; xlabel="Zeit [s]", ylabel="|V₁| [pu]", title="Spannungsbetrag Bus A (Knoten 1): |V₁| = √(u_r²+u_i²)")
        lines!(ax1, t_ref, ref_fullsim.bus1_V[idx];      label="PF  bus1_V", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax1, ts,    sqrt.(ur_jl.^2 .+ ui_jl.^2); label="Julia  √(u_r²+u_i²)", color=:red,  linewidth=2, linestyle=:dash)
        axislegend(ax1; position=:rb)

        ax2 = Axis(fig[1,2]; xlabel="Zeit [s]", ylabel="δ₁ [°]", title="Spannungswinkel Bus A (Knoten 1): δ₁ = atan(u_i, u_r)")
        lines!(ax2, t_ref, ref_fullsim.bus1_delta[idx];     label="PF  bus1_delta", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax2, ts,    atan.(ui_jl, ur_jl) .* (180/π); label="Julia  atan(u_i, u_r) [°]", color=:red,  linewidth=2, linestyle=:dash)
        axislegend(ax2; position=:rb)

        ax3 = Axis(fig[2,1]; xlabel="Zeit [s]", ylabel="u_r [pu]", title="Realteil der Spannung Bus A (Knoten 1): u_r = Re{V₁}")
        lines!(ax3, t_ref, ref_fullsim.bus1_vr[idx]; label="PF  bus1_vr", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax3, ts,    ur_jl;                    label="Julia  busbar.u_r", color=:red,  linewidth=2, linestyle=:dash)
        axislegend(ax3; position=:rb)

        ax4 = Axis(fig[2,2]; xlabel="Zeit [s]", ylabel="u_i [pu]", title="Imaginärteil der Spannung Bus A (Knoten 1): u_i = Im{V₁}")
        lines!(ax4, t_ref, ref_fullsim.bus1_vi[idx]; label="PF  bus1_vi", color=:blue, linewidth=2, alpha=0.7)
        lines!(ax4, ts,    ui_jl;                    label="Julia  busbar.u_i", color=:red,  linewidth=2, linestyle=:dash)
        axislegend(ax4; position=:rb)

        fig
    end
    save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","PowerFactory_valid","PV_3bus_fullsim_busA.png"), fig_fullsim_busA)
end

# Voltage-Dip-Flags und Limiter-Zustände (nur Julia-intern, kein PF-Vergleich)
if isdefined(Main, :EXPORT_FIGURES) && Main.EXPORT_FIGURES
    fig_fullsim_vdip = let
        ts = let r = refine_timeseries(sol_fullsim.t, 20); r[(r .>= 0) .& (r .<= 0.3)] end
        fig = Figure(resolution=(700, 400))
        ax  = Axis(fig[1,1]; xlabel="Zeit [s]", ylabel="[-]", title="Voltage-Dip-Flags und Limiter-Zustände (Julia-intern, kein PF-Vergleich)")
        lines!(ax, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊reecb₊Voltage_dip)).u; label="reecb.Voltage_dip  (REEC-B Spannungseinbruch-Flag)", color=:blue,      linewidth=2)
        lines!(ax, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊repca₊Voltage_dip)).u; label="repca.Voltage_dip  (REPC-A Spannungseinbruch-Flag)", color=:red,       linewidth=2)
        lines!(ax, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊repca₊Q_lim)).u;       label="repca.Q_lim         (REPC-A Q-Begrenzungs-Flag)",    color=:green,     linewidth=2)
        lines!(ax, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊reecb₊P_lim)).u;       label="reecb.P_lim         (REEC-B P-Begrenzungs-Flag)",    color=:pink,      linewidth=2)
        lines!(ax, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊regca₊I_qr)).u;        label="regca.I_qr          (REGC-A Blindstrom-Rampenrate)", color=:turquoise, linewidth=2)
        lines!(ax, ts, sol_fullsim(ts, idxs=VIndex(:GEN1, :PV₊regca₊I_pr)).u;        label="regca.I_pr          (REGC-A Wirkstrom-Rampenrate)",  color=:orange,    linewidth=2)
        fig
    end
    save(joinpath(pkgdir(OpPoDyn),"docs","src","assets","PowerFactory_valid","PV_3bus_fullsim_vdip.png"), fig_fullsim_vdip)
end

# Konsistenzcheck: erste Werte nach t = 0.1 s an Bus A
let
    println("\n", "="^70)
    println("Erste Werte nach t = 0.1 s — Bus A (vollständige Simulation)")
    println("="^70)
    csv_idx = findfirst(ref_fullsim.time .> 0.1)
    println("\nPowerFactory (ref_fullsim) — t = $(ref_fullsim.time[csv_idx]) s")
    println("  |V| = $(ref_fullsim.bus1_V[csv_idx])  δ = $(ref_fullsim.bus1_delta[csv_idx]) °")
    println("  u_r = $(ref_fullsim.bus1_vr[csv_idx])  u_i = $(ref_fullsim.bus1_vi[csv_idx])")
    jl_t   = sol_fullsim.t[findfirst(sol_fullsim.t .> 0.1)]
    ur_jl1 = sol_fullsim(jl_t, idxs=VIndex(1, :busbar₊u_r))
    ui_jl1 = sol_fullsim(jl_t, idxs=VIndex(1, :busbar₊u_i))
    println("\nJulia (sol_fullsim)        — t = $jl_t s")
    println("  |V| = $(sqrt(ur_jl1^2+ui_jl1^2))  δ = $(atan(ui_jl1, ur_jl1)*180/π) °")
    println("  u_r = $ur_jl1  u_i = $ui_jl1")
end


# ══════════════════════════════════════════════════════════════════════════════
#  Open-Loop-Test mit vorgegebenen Klemmengrößen (Prescribed Inputs)
#   Spannungen und Ströme werden als Inputs aus PF-CSV vorgegeben.
#   REGCA/REEC/REPC laufen frei; pir_out / pii_out sind Modellausgang.
# ══════════════════════════════════════════════════════════════════════════════

include(joinpath(@__DIR__, "open_loop_utils.jl"))
include(joinpath(@__DIR__, "prescribed_pv_test.jl"))

# Gemeinsame Diagnose-Schlüssel für comparison_figure
regc_keys        = ["regc_Iqr", "regc_o2", "regc_o3", "regc_dIq", "regc_Qgen", "regc_Vtfiltlim"]
repc_q_diag_keys = ["repc_Qfltr", "repc_dQ", "repc_dQin", "repc_dQdbd", "repc_Qe", "repc_Qlim"]
exclude_keys     = union(regc_keys, repc_q_diag_keys)
repc_q_keys      = ["repc_Qbranch", "repc_Qfltr", "repc_dQ", "repc_dQin",
                    "repc_dQdbd", "repc_Qe", "repc_Qlim", "repc_Qext"]


# ── 2a: Mit Switch-off-Schwellwert (PF-Standard: 0.1 pu / 0.15 pu) ───────────

# Referenz Klemmengrößen: PF mit Switch-off-Schwellwert 0.1 pu / 0.15 pu
ref_inputs_thresh = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV_pf",
        "testcase3Bus-with-event_InputMeasures_moretimesteps.csv"),
    DataFrame;
    header = 3,
    delim = ';',
    decimal = ',',
    silencewarnings = true
)
rename!(ref_inputs_thresh, 1 => :time)

# Referenz interne Variablen: PF mit Switch-off-Schwellwert 0.1 pu / 0.15 pu
ref_vars_thresh = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV_pf",
        "variables-testcase-with-event_moretimesteps.csv"),
    DataFrame;
    header = 3,
    decimal = ',',
    drop = (i, name) -> contains(string(name), "nrows="),
    silencewarnings = true
)

# Lösung: prescribed inputs, Referenz mit Switch-off-Schwellwert 0.1/0.15 pu
sol_prescribed_thresh = prescribed_pv_test(
    ref_inputs_thresh;
    obs_df = ref_vars_thresh,
)

# Hauptvergleich (ohne REGC-Interna und REPC-Q-Diagnose)
fig_prescribed_thresh = comparison_figure(
    (; sol=sol_prescribed_thresh.sol,
       dat=OrderedDict(k => v for (k, v) in sol_prescribed_thresh.dat if k ∉ exclude_keys));
    tmin=0.0, tmax=0.3,
)

# REGC-Interna: I_qr, o2, o3, ΔI_q, Q_gen, V_tfiltlim
fig_prescribed_thresh_regc = comparison_figure(
    (; sol=sol_prescribed_thresh.sol,
       dat=OrderedDict(k => sol_prescribed_thresh.dat[k] for k in regc_keys));
    tmin=0.0, tmax=0.3,
)

# REPC Q-Signalpfad: Q_branch → Q_fltr → ΔQ → ΔQ_in → ΔQ_dbd → Q_e → Q_lim → Q_ext
fig_prescribed_thresh_repc_q = comparison_figure(
    (; sol=sol_prescribed_thresh.sol,
       dat=OrderedDict(k => sol_prescribed_thresh.dat[k] for k in repc_q_keys));
    tmin=0.0, tmax=0.3,
)

# CSV-Konsistenzcheck: √(pir²+pii²) vs repc_Ibranch, pvi·pir−pvr·pii vs repc_Qbranch
fig_prescribed_thresh_csvcheck = let
    idx  = (ref_inputs_thresh.time .>= 0.09) .& (ref_inputs_thresh.time .<= 0.14)
    t    = ref_inputs_thresh.time[idx]
    pir  = ref_inputs_thresh.pir[idx];  pii = ref_inputs_thresh.pii[idx]
    pvr  = ref_inputs_thresh.pvr[idx];  pvi = ref_inputs_thresh.pvi[idx]
    idx2 = (ref_vars_thresh.time .>= 0.09) .& (ref_vars_thresh.time .<= 0.14)
    t2   = ref_vars_thresh.time[idx2]

    fig = Figure(size=(900, 400))
    ax1 = Axis(fig[1,1]; xlabel="Time [s]", title="I_branch: √(pir²+pii²) vs repc_Ibranch")
    lines!(ax1, t,  sqrt.(pir.^2 .+ pii.^2);           label="√(pir²+pii²) [ref_inputs]",  linewidth=2)
    lines!(ax1, t2, ref_vars_thresh.repc_Ibranch[idx2]; label="repc_Ibranch [ref_vars]",     linewidth=2, linestyle=:dash)
    lines!(ax1, t,  ref_inputs_thresh.I_measure[idx];   label="I_measure    [ref_inputs]",   linewidth=2, color=:black)
    axislegend(ax1; position=:rt)

    ax2 = Axis(fig[1,2]; xlabel="Time [s]", title="Q_branch: pvi·pir−pvr·pii vs repc_Qbranch")
    lines!(ax2, t,  pvi.*pir .- pvr.*pii;               label="pvi·pir−pvr·pii [ref_inputs]", linewidth=2)
    lines!(ax2, t2, ref_vars_thresh.repc_Qbranch[idx2]; label="repc_Qbranch    [ref_vars]",   linewidth=2, linestyle=:dash)
    lines!(ax2, t,  ref_inputs_thresh.Q_measure[idx];   label="Q_measure       [ref_inputs]", linewidth=2, color=:black)
    axislegend(ax2; position=:rt)
    fig
end

# Q_fltr Detailvergleich: Q_branch, simpleLag2-Zustand, Q_fltr
# Referenz: separate hochauflösende PF-CSV mit Switch-off-Schwellwert
ref_qfltr_thresh = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV_pf",
        "testcase3Bus-with-event_QbranchQfltr_moretimesteps.csv"),
    DataFrame;
    header = 3,
    delim = ';',
    decimal = ',',
    silencewarnings = true
)
rename!(ref_qfltr_thresh, 1 => :time)

fig_prescribed_thresh_qfltr = let
    sys   = sol_prescribed_thresh.sys
    sol   = sol_prescribed_thresh.sol
    tsim  = _ol_refine(sol.t, 5)
    tsim  = tsim[(tsim .>= 0.0) .& (tsim .<= 0.3)]
    idx   = (ref_qfltr_thresh.time .>= 0.0) .& (ref_qfltr_thresh.time .<= 0.3)
    t_ref = ref_qfltr_thresh.time[idx]
    idx_in = (ref_inputs_thresh.time .>= 0.0) .& (ref_inputs_thresh.time .<= 0.3)

    fig = Figure(size=(900, 900))
    ax1 = Axis(fig[1,1]; xlabel="Zeit [s]", ylabel="Q [pu]", title="Blindleistung Q am Einspeisepunkt (Eingang REPC-A Q-Filter)")
    lines!(ax1, ref_inputs_thresh.time[idx_in], ref_inputs_thresh.Q_measure[idx_in]; label="PF  Q_measure (ref_inputs_thresh)", color=:blue, linewidth=2, alpha=0.7)
    lines!(ax1, tsim, sol(tsim; idxs=sys.Q_measure).u; label="Julia  Q_measure", color=:red, linewidth=2, linestyle=:dash)
    axislegend(ax1)

    ax2 = Axis(fig[2,1]; xlabel="Zeit [s]", ylabel="[pu]", title="Filterzustand x_fltr2 (simpleLag2-Ausgang = Q_fltr, Eingang: Q_branch)")
    lines!(ax2, t_ref, ref_qfltr_thresh.xfltr2[idx]; label="PF  xfltr2 (ref_qfltr_thresh)", color=:blue, linewidth=2, alpha=0.7)
    lines!(ax2, tsim, sol(tsim; idxs=sys.repca.simpleLag2.out).u; label="Julia  repca.simpleLag2.out (= Q_fltr)", color=:red, linewidth=2, linestyle=:dash)
    axislegend(ax2)

    ax3 = Axis(fig[3,1]; xlabel="Zeit [s]", ylabel="Q_fltr [pu]", title="Gefilterter Blindleistungswert Q_fltr (Ausgang REPC-A Q-Filter)")
    lines!(ax3, t_ref, ref_qfltr_thresh.Q_fltr[idx]; label="PF  Q_fltr (ref_qfltr_thresh)", color=:blue, linewidth=2, alpha=0.7)
    lines!(ax3, tsim, sol(tsim; idxs=sys.repca.Q_fltr).u; label="Julia  repca.Q_fltr", color=:red, linewidth=2, linestyle=:dash)
    axislegend(ax3)
    fig
end

# ── 2b: Mit berechneten Measures (aus pvr/pvi/pir/pii statt CSV) ─────────────

sol_prescribed_thresh_calc = prescribed_pv_test(
    ref_inputs_thresh;
    obs_df        = ref_vars_thresh,
    calc_measures = true,
)

fig_prescribed_thresh_calc = comparison_figure(
    (; sol=sol_prescribed_thresh_calc.sol,
       dat=OrderedDict(k => v for (k, v) in sol_prescribed_thresh_calc.dat if k ∉ exclude_keys));
    tmin=0.0, tmax=0.3,
)

fig_prescribed_thresh_calc_repc_q = comparison_figure(
    (; sol=sol_prescribed_thresh_calc.sol,
       dat=OrderedDict(k => sol_prescribed_thresh_calc.dat[k] for k in repc_q_keys));
    tmin=0.0, tmax=0.3,
)

# ──  Ohne Switch-off-Schwellwert (PF-Schwellwert = 0 pu, deaktiviert) ──────

# Referenz Klemmengrößen: PF ohne Switch-off-Schwellwert (0 pu, deaktiviert)
ref_inputs_nothresh = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV_pf",
        "testcase3Bus-with-event_InputMeasures_withoutThresholdOperationVoltage.csv"),
    DataFrame;
    header = 3,
    delim = ';',
    decimal = ',',
    silencewarnings = true
)
rename!(ref_inputs_nothresh, 1 => :time)

# Referenz interne Variablen: PF ohne Switch-off-Schwellwert (0 pu, deaktiviert)
ref_vars_nothresh = CSV.read(
    joinpath(pkgdir(OpPoDyn),"test","WECC_model_tests","PV_pf",
        "variables-testcase3Bus-with-event_addVariableNames_withoutThresholdOperationVoltage.csv"),
    DataFrame;
    header = 3,
    decimal = ',',
    drop = (i, name) -> contains(string(name), "nrows="),
    silencewarnings = true
)

# Lösung: prescribed inputs, Referenz ohne Switch-off-Schwellwert
sol_prescribed_nothresh = prescribed_pv_test(
    ref_inputs_nothresh;
    obs_df = ref_vars_nothresh,
)

# Hauptvergleich (ohne REGC-Interna und REPC-Q-Diagnose)
fig_prescribed_nothresh = comparison_figure(
    (; sol=sol_prescribed_nothresh.sol,
       dat=OrderedDict(k => v for (k, v) in sol_prescribed_nothresh.dat if k ∉ exclude_keys));
    tmin=0.0, tmax=0.3,
)

# REGC-Interna: I_qr, o2, o3, ΔI_q, Q_gen, V_tfiltlim
fig_prescribed_nothresh_regc = comparison_figure(
    (; sol=sol_prescribed_nothresh.sol,
       dat=OrderedDict(k => sol_prescribed_nothresh.dat[k] for k in regc_keys));
    tmin=0.0, tmax=0.3,
)

# REPC Q-Signalpfad: Q_branch → Q_fltr → ΔQ → ΔQ_in → ΔQ_dbd → Q_e → Q_lim → Q_ext
fig_prescribed_nothresh_repc_q = comparison_figure(
    (; sol=sol_prescribed_nothresh.sol,
       dat=OrderedDict(k => sol_prescribed_nothresh.dat[k] for k in repc_q_keys));
    tmin=0.0, tmax=0.3,
)

# ══════════════════════════════════════════════════════════════════════════════
# Ohne Schwellwert, Measures berechnet (calc_measures=true):
# I/P/Q_measure werden aus pvr/pvi/pir/pii berechnet, nicht direkt aus CSV.
# Vergleich gegen PF-Referenz ohne Schwellwert.
# ══════════════════════════════════════════════════════════════════════════════

sol_prescribed_nothresh_calc = prescribed_pv_test(
    ref_inputs_nothresh;
    obs_df        = ref_vars_nothresh,
    calc_measures = true,
)

# ── Plot 1: Inputs + berechnete Measures vs PF-Referenz ──────────────────────
fig_nothresh_calc_inputs = let
    sys  = sol_prescribed_nothresh_calc.sys
    sol  = sol_prescribed_nothresh_calc.sol
    tsim = _ol_refine(sol.t, 5)
    tsim = tsim[(tsim .>= 0.0) .& (tsim .<= 0.3)]

    idx_in = (ref_inputs_nothresh.time .>= 0.0) .& (ref_inputs_nothresh.time .<= 0.3)
    t_in   = ref_inputs_nothresh.time[idx_in]

    Vt   = sol(tsim; idxs=sys.V_t).u
    δv   = sol(tsim; idxs=sys.δ_v).u

    fig = Figure(size=(1200, 1100))

    ax1 = Axis(fig[1,1]; xlabel="Zeit [s]", ylabel="[pu]", title="pvr (Spannung real)")
    lines!(ax1, t_in,  ref_inputs_nothresh.pvr[idx_in]; label="PF pvr",          color=:blue, linewidth=2, alpha=0.7)
    lines!(ax1, tsim,  Vt .* cos.(δv);                  label="Julia V_t·cos(δ)", color=:red,  linewidth=2, linestyle=:dash)
    axislegend(ax1)

    ax2 = Axis(fig[1,2]; xlabel="Zeit [s]", ylabel="[pu]", title="pvi (Spannung imag)")
    lines!(ax2, t_in,  ref_inputs_nothresh.pvi[idx_in]; label="PF pvi",          color=:blue, linewidth=2, alpha=0.7)
    lines!(ax2, tsim,  Vt .* sin.(δv);                  label="Julia V_t·sin(δ)", color=:red,  linewidth=2, linestyle=:dash)
    axislegend(ax2)

    ax3 = Axis(fig[2,1]; xlabel="Zeit [s]", ylabel="[pu]", title="pir (Strom real, Eingabe)")
    lines!(ax3, t_in,  ref_inputs_nothresh.pir[idx_in]; label="PF pir (Eingabe)", color=:blue, linewidth=2)
    axislegend(ax3)

    ax4 = Axis(fig[2,2]; xlabel="Zeit [s]", ylabel="[pu]", title="pii (Strom imag, Eingabe)")
    lines!(ax4, t_in,  ref_inputs_nothresh.pii[idx_in]; label="PF pii (Eingabe)", color=:blue, linewidth=2)
    axislegend(ax4)

    ax5 = Axis(fig[3,1]; xlabel="Zeit [s]", ylabel="[pu]", title="I_measure: Julia (berechnet) vs PF")
    lines!(ax5, t_in,  ref_inputs_nothresh.I_measure[idx_in]; label="PF I_measure (CSV)",       color=:blue,  linewidth=2, alpha=0.7)
    lines!(ax5, tsim,  sol(tsim; idxs=sys.I_measure).u;       label="Julia √(pir²+pii²)",        color=:red,   linewidth=2, linestyle=:dash)
    axislegend(ax5)

    ax6 = Axis(fig[3,2]; xlabel="Zeit [s]", ylabel="[pu]", title="P_measure: Julia (berechnet) vs PF")
    lines!(ax6, t_in,  ref_inputs_nothresh.P_measure[idx_in]; label="PF P_measure (CSV)",        color=:blue,  linewidth=2, alpha=0.7)
    lines!(ax6, tsim,  sol(tsim; idxs=sys.P_measure).u;       label="Julia pvr·pir+pvi·pii",      color=:red,   linewidth=2, linestyle=:dash)
    axislegend(ax6)

    ax7 = Axis(fig[4,1]; xlabel="Zeit [s]", ylabel="[pu]", title="Q_measure: Julia (berechnet) vs PF")
    lines!(ax7, t_in,  ref_inputs_nothresh.Q_measure[idx_in]; label="PF Q_measure (CSV)",        color=:blue,  linewidth=2, alpha=0.7)
    lines!(ax7, tsim,  sol(tsim; idxs=sys.Q_measure).u;       label="Julia pvi·pir−pvr·pii",      color=:red,   linewidth=2, linestyle=:dash)
    axislegend(ax7)

    fig
end

# ── Plot 2: Modulergebnisse (REGC/REEC/REPC) vs PF-Referenz ohne Schwellwert ─
fig_nothresh_calc_modules = comparison_figure(
    (; sol=sol_prescribed_nothresh_calc.sol,
       dat=OrderedDict(k => v for (k, v) in sol_prescribed_nothresh_calc.dat if k ∉ exclude_keys));
    tmin=0.0, tmax=0.3,
)

fig_nothresh_calc_repc_q = comparison_figure(
    (; sol=sol_prescribed_nothresh_calc.sol,
       dat=OrderedDict(k => sol_prescribed_nothresh_calc.dat[k] for k in repc_q_keys));
    tmin=0.0, tmax=0.3,
)

# ── CSV-Konsistenzcheck: PF-Measures vs. aus pvr/pvi/pir/pii berechnete Werte ─
# Reine Datenvergleich, keine Simulation — zeigt ob PF intern filtert/abweicht.
fig_csv_measures_check = let
    idx = (ref_inputs_nothresh.time .>= 0.0) .& (ref_inputs_nothresh.time .<= 0.3)
    t   = ref_inputs_nothresh.time[idx]
    pvr = ref_inputs_nothresh.pvr[idx]
    pvi = ref_inputs_nothresh.pvi[idx]
    pir = ref_inputs_nothresh.pir[idx]
    pii = ref_inputs_nothresh.pii[idx]

    I_calc = sqrt.(pir.^2 .+ pii.^2)
    P_calc = pvr.*pir .+ pvi.*pii
    Q_calc = pvi.*pir .- pvr.*pii

    fig = Figure(size=(900, 1200))

    ax1 = Axis(fig[1,1]; xlabel="Zeit [s]", ylabel="[pu]", title="pir (Strom real)")
    lines!(ax1, t, pir; label="PF pir (CSV)", color=:blue, linewidth=2)
    axislegend(ax1)

    ax2 = Axis(fig[2,1]; xlabel="Zeit [s]", ylabel="[pu]", title="pii (Strom imag)")
    lines!(ax2, t, pii; label="PF pii (CSV)", color=:blue, linewidth=2)
    axislegend(ax2)

    ax3 = Axis(fig[3,1]; xlabel="Zeit [s]", ylabel="[pu]", title="I_measure: CSV vs √(pir²+pii²)")
    lines!(ax3, t, ref_inputs_nothresh.I_measure[idx]; label="PF I_measure (CSV)", color=:blue, linewidth=2)
    lines!(ax3, t, I_calc;                             label="√(pir²+pii²)",       color=:red,  linewidth=2, linestyle=:dash)
    axislegend(ax3)

    ax4 = Axis(fig[4,1]; xlabel="Zeit [s]", ylabel="[pu]", title="P_measure: CSV vs pvr·pir+pvi·pii")
    lines!(ax4, t, ref_inputs_nothresh.P_measure[idx]; label="PF P_measure (CSV)", color=:blue, linewidth=2)
    lines!(ax4, t, P_calc;                             label="pvr·pir+pvi·pii",     color=:red,  linewidth=2, linestyle=:dash)
    axislegend(ax4)

    ax5 = Axis(fig[5,1]; xlabel="Zeit [s]", ylabel="[pu]", title="Q_measure: CSV vs pvi·pir−pvr·pii")
    lines!(ax5, t, ref_inputs_nothresh.Q_measure[idx]; label="PF Q_measure (CSV)", color=:blue, linewidth=2)
    lines!(ax5, t, Q_calc;                             label="pvi·pir−pvr·pii",     color=:red,  linewidth=2, linestyle=:dash)
    axislegend(ax5)

    fig
end
