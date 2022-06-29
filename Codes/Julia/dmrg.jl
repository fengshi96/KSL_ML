using SparseArrays
using LinearAlgebra
using ITensors
using HDF5

include("src/Lattice.jl")
include("src/Parameter.jl")


para = GetParameter("input.inp")
Nsite, mesh, nn, _, _ = Honeycomb(para)
Nstates = para.Nstates
Kx = para.Kxx
Ky = para.Kyy
Kz = para.Kzz
Hx = para.Hx
Hy = para.Hy
Hz = para.Hz

# @time Ham = Hamiltonian(param)
# @time eigvals, eigvecs = eigs(Ham, nev = Nstates, which=:SR, tol = 0.001) # SR = Smallest Real 
# println("\nEigen values:"); show(stdout, "text/plain", real.(eigvals)); println()

let
    sites = siteinds("S=1/2", Nsite)
    ampo = OpSum()

    for j=1:Nsite
        nx::Int64 = nn[j, 1]
        ny::Int64 = nn[j, 2]
        nz::Int64 = nn[j, 3]
        
        if j < nn[j, 1]
            ampo += Kx, "Sx", j, "Sx", nx
        end
        if j < nn[j, 2]
            ampo += Ky, "Sy", j, "Sy", ny
        end
        if j < nn[j, 3]
            ampo += Kz, "Sz", j, "Sz", nz
        end
    end
    
    H = MPO(ampo, sites)
    psi0 = randomMPS(sites, 100)

    sweeps = Sweeps([
         "maxdim" "mindim" "cutoff" "noise"
          100      100       1e-12    1E-7
          200      100       1e-12    1E-8
          200      100       1e-12    0
          200      100       1e-12    0
          300      100       1e-12    0
          300      100       1e-12    0
          300      100       1e-12    0
          300      100       1e-12    0
         ])

    energy, psi = dmrg(H, psi0, sweeps)

    f = h5open("data.h5","w")
    write(f,"psi",psi)
    close(f)

end





