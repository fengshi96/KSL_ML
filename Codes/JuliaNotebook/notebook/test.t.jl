using Random
using LinearAlgebra
using ITensors 
using HDF5 

f = h5open("data.h5","r")
psi = read(f,"psi",MPS)
close(f)
