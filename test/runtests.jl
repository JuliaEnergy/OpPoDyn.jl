using Test
using SafeTestsets
using Aqua
using ExplicitImports
using NetworkDynamics
using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using Makie
using CairoMakie
using OrderedCollections

using OpPoDyn
using OpPoDyn.Library
using PowerDynamics: PowerDynamics

@testset "OpPoDyn.jl Tests" begin
    @testset "Package Quality Tests" begin
        @info "Begin Package quality tests"
        Aqua.test_all(OpPoDyn;
            ambiguities=false,
            persistent_tasks=false)
        @test_broken isempty(Docs.undocumented_names(OpPoDyn))

        @test check_no_implicit_imports(OpPoDyn; skip=(Base, Core, NetworkDynamics, PowerDynamics)) === nothing
        # mtkmodel macro depends on some symbols
        @test_broken check_no_stale_explicit_imports(OpPoDyn) === nothing

        path = joinpath(pkgdir(OpPoDyn),"src","Library","Library.jl")
        @test check_no_implicit_imports(OpPoDyn.Library, path) === nothing
        # mtkmodel macro depends on some symbols
        @test_broken check_no_stale_explicit_imports(OpPoDyn.Library, path) === nothing
    end

    @safetestset "Library tests" begin include("Library_test.jl") end
    @safetestset "utils tests" begin include("utils_test.jl") end
end
