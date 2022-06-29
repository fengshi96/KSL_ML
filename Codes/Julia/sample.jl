using LinearAlgebra
using ITensors
using HDF5
using Random
include("src/sampler.jl")

let
    Nsite = 56
    f = h5open("data/data_TC.h5","r")
    @time ψ = read(f,"psi",MPS)
    close(f)
    siteSet = siteinds(ψ)   
#     siteSet = siteinds("S=1/2", Nsite)
#     ψ = replace_siteinds(ψ,siteSet)

    sites = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    axes = ["z", "z", "z", "z", "z", "z", "z", "z", "z", "z", "z", "z", "z", "z", "z", "z"]
    @time sample, probs, clamps = snapShot(sites, ψ, axes)#, rand(100:9000))
    print(clamps)

    data = hcat(sites, axes, sample, probs)
    println(); show(stdout, "text/plain", data); println()

#     orthogonalize!(ψ, 2)
#     println(ψ[2])
end
