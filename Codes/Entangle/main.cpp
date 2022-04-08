#include <iostream>
#include "itensor/all.h"
#include <Eigen/Dense>
#include <iomanip>

using namespace itensor;

int main(int argc, char* argv[]) {
    if(argc < 2)  {
        printfln("Usage: %s input_file",argv[0]);
        return 0;
    }
    auto input = InputGroup(argv[1],"input");
    int Nx = input.getInt("Nx");
    int Ny = input.getInt("Ny");
    int N = Nx*Ny*2;
    std::string sysIndx = input.getString("sysIndx");
    int cutLabel = 1;
    for (char const &c : sysIndx) if(c == ',') cutLabel++;
    std::cout << "cutLabel = " << cutLabel << std::endl;


    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    auto ampo = AutoMPO(sites);
    std::cout << " Reading from File " << std::endl;
    readFromFile("sites.dat",sites);
    auto psi = readFromFile<MPS>(std::string("psi.dat"),sites);


    // Begin entanglment
    int bond = cutLabel;  // index of which bond to cut
    psi.position(bond);

    //SVD this wavefunction to get the spectrum of density-matrix eigenvalues
    auto l = leftLinkIndex(psi,bond);
    auto s = siteIndex(psi,bond);
    auto [U,S,V] = svd(psi(bond),{l,s});
    auto u = commonIndex(U,S);

    //Apply von Neumann formula
    //to the squares of the singular values
    Real SvN = 0.;
    for(auto n : range1(dim(u)))
    {
        auto Sn = elt(S,n,n);
        auto p = sqr(Sn);
        if(p > 1E-12) SvN += -p*log(p);
    }
    printfln("Across bond b=%d, SvN = %.10f",bond,SvN);


    return 0;


}