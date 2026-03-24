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


v1a = VertexModel(pfPQ(name = :bus1, P = 0.88, Q = -0.33; current_source=true), vidx=1)
v2 = VertexModel(pfPV(name= :junction, P=0, V=1), vidx=2)
v3 = VertexModel(pfSlack(name = :slack, V = 1, δ=0), vidx=3)
e1 = LoopbackConnection(; src=:bus1, dst=:junction, potential=[:u_r, :u_i], flow=[:i_r, :i_i])
e2 =  compile_line(MTKLine(PiLine(; name=:PwLine)),src=:junction, dst=:slack)

pfnw = Network([v1a, v2, v3], [e1, e2])
solve_powerflow(pfnw) #works



v1b = VertexModel(pfPQ(name = :bus1, P = 0.88, Q = -0.33; current_source=true, assume_io_coupling=true), vidx=1)
pfnw_currentsource = Network([v1b, v2, v3], [e12, e23])
solve_powerflow(pfnw_currentsource) #Error