using IntervalCensoredMultistate
using Test

@testset "IntervalCensoredMultistate.jl" begin
    exported = names(IntervalCensoredMultistate)
    @test :fit_single_event in exported
    @test :fit_competing_risks in exported
    @test :fit_multistate in exported
    @test :reg_est_bone in exported
end
