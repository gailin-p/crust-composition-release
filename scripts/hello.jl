# examples/01-hello.jl
using MPI
MPI.Init()

comm = MPI.COMM_WORLD
print("Hello world, I am rank $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))\n")
MPI.Barrier(comm)