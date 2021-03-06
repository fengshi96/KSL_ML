#include <iostream>
#include "itensor/all.h"
#include "src/honeycomb.h"
#include "src/operators.h"

using namespace itensor;

int main(int argc, char* argv[]) {
    if(argc < 2)  {
        printfln("Usage: %s input_file",argv[0]);
        return 0;
    }
    auto input = InputGroup(argv[1],"input");
    int Nx = input.getInt("Nx");
    int Ny = input.getInt("Ny");
    std::string Coupling = input.getString("Coupling");
    int N = Nx*Ny*2;
    bool xperiodic = input.getInt("IsPeriodicX");
    bool yperiodic = input.getInt("IsPeriodicY");

    double Kx = input.getReal("Kx");
    double Ky = input.getReal("Ky");
    double Kz = input.getReal("Kz");

    double Hx = input.getReal("Hx");
    double Hy = input.getReal("Hy");
    double Hz = input.getReal("Hz");

    // sweep parameters
    auto totalSweeps = input.getInt("totalSweeps");
    auto sw_table = InputGroup(input,"sweep_table");
    auto sweeps = Sweeps(totalSweeps,sw_table);
    println(sweeps);


    //  Initialize the site degrees of freedom.
    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    // auto lattice = honeycombLattice(Nx,Ny,{"YPeriodic=",yperiodic, "XPeriodic=",xperiodic, "Cutting=",cutting});
    auto lattice = honeycombLattice(Nx,Ny,{"YPeriodic=",yperiodic, "XPeriodic=",xperiodic});
    std::cout << lattice;

    auto ampo = AutoMPO(sites);
    if (Coupling == "Kitaev") {
        // Kitaev interaction
        for(auto bnd : lattice) {
            if (bnd.type == "Sx")
                ampo +=  Kx, bnd.type, bnd.s1, bnd.type, bnd.s2;
            else if (bnd.type == "Sy")
                ampo += Ky, bnd.type, bnd.s1, bnd.type, bnd.s2;
            else if (bnd.type == "Sz")
                ampo += Kz, bnd.type, bnd.s1, bnd.type, bnd.s2;
        }
    } else if (Coupling == "Heisenberg") {
        // Heisenberg interaction
        for(auto bnd : lattice) {
            ampo +=  Kx, "Sx", bnd.s1, "Sx", bnd.s2;
            ampo +=  Ky, "Sy", bnd.s1, "Sy", bnd.s2;
            ampo +=  Kz, "Sz", bnd.s1, "Sz", bnd.s2;
        }
    }



    // Add magnetic field
    for (int i = 1; i <=N; ++i) {
        ampo += Hx, "Sx", i;
        ampo += Hy, "Sy", i;
        ampo += Hz, "Sz", i;
    }

    auto H = toMPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i)
    {
        if(i%2 == 1)
            state.set(i,"Up");
        else
            state.set(i,"Dn");
    }

    auto psi0 = MPS(state);

    // overlap calculates matrix elements of MPO's with respect to MPS's
    // inner(psi0,H,psi0) = <psi0|H|psi0>
    printfln("Initial energy = %.5f", innerC(psi0,H,psi0) );

    // Begin the DMRG calculation
    auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing overlap = %.10f", innerC(psi,H,psi) );
    writeToFile(std::string("sites.dat"),sites);
    writeToFile(std::string("psi.dat"),psi);

    //auto f = h5_open()


    return 0;
}
