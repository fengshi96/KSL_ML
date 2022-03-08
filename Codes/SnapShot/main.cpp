#include <iostream>
#include "itensor/all.h"
#include <Eigen/Dense>
#include <iomanip>
#include <array>
#include <algorithm>

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

    // ================= Read in mask ====================================
    int maskSize = input.getInt("maskSize");
    std::vector<int> mpsIndx(maskSize);
    std::vector<int> x(maskSize);
    std::vector<int> y(maskSize);
    std::vector<int> z(maskSize);
    auto mask = InputGroup(input, "mask");
    printfln("mpsId x y z");
    if (mask.GotoGroup()) {
        mask.SkipLine();
        for (int i = 0; i<maskSize; i++) {
            mask.file() >> mpsIndx[i] >> x[i] >> y[i] >> z[i];
            printfln("  %d   %d %d %d",mpsIndx[i],x[i],y[i],z[i]);
        }
    } else {
        Error("Couldn't find mask " + mask.name());
    }


    // ============== Read wave function and initialize MPS =================
    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    auto ampo = AutoMPO(sites);
    std::cout << " Reading from File " << std::endl;
    readFromFile("sites.dat",sites);
    auto psi = readFromFile<MPS>(std::string("psi.dat"),sites);

    // =================== Multi-Point Observables =========================
    int start = *std::min_element(mpsIndx.begin(), mpsIndx.end()); println("Mask start at site: ", start);
    int end = *std::max_element(mpsIndx.begin(), mpsIndx.end()); println("Mask ends at site: ", end);

    // initialize the tensor train from the left side
    psi.position(start);
    ITensor C = psi(start);

    // the head of tensor train
    std::string S_op = (x[start] == 1 ? "Sx" : (y[start] == 1 ? "Sy" : "Sz"));
    C *= op(sites, S_op, start);
    auto ir = commonIndex(psi(start),psi(start+1),"Link");
    C *= dag(prime(prime(psi(start), "Site"), ir));
    C *= psi(start+1);

    // the bulk of tensor train
    for (int i = start+1; i < end; i++) { // i-th neighbor to the right of the head
        if (std::find(mpsIndx.begin(), mpsIndx.end(), i) != mpsIndx.end()) {
            // if i+1 is in mpsIndx (continuous indices)
            S_op = (x[i] == 1 ? "Sx" : (y[i] == 1 ? "Sy" : "Sz"));
            C *= op(sites, S_op, i);
            C *= dag(prime(prime(psi(i), "Site"), "Link"));
            C *= psi(i+1);
        } else {
            C *= dag(prime(psi(i), "Link"));
            C *= psi(i+1);
        }
    }

    // the tail of tensor train
    S_op = (x[end] == 1 ? "Sx" : (y[end] == 1 ? "Sy" : "Sz"));
    C *= op(sites, S_op, end);
    auto il = commonIndex(psi(end-1),psi(end),"Link");
    C *= dag(prime(prime(psi(end), "Site"), il));

    auto results = eltC(C);
    std::cout << "Result = " << results << std::endl;

    return 0;
}

//    // ====test====
//    ITensor Sz_1 = op(sites, "Sz", 1);
//    ITensor Sz_2 = op(sites, "Sz", 2);
//    ITensor Sz_3 = op(sites, "Sz", 3);
//
//    psi.position(1);
//
//    ITensor C = psi(1);
//    C *= Sz_1;  // with primed physical index
//    auto ir = commonIndex(psi(1),psi(2),"Link");
//    C *= dag(prime(prime(psi(1), "Site"), ir));
//
//    C *= psi(2);
//    C *= Sz_2;
//    C *= dag(prime(prime(psi(2), "Site"), "Link"));
//
//    C *= psi(3);
//    C *= Sz_3;
//    auto il = commonIndex(psi(2),psi(3),"Link");
//    C *= dag(prime(prime(psi(3), "Site"), il));
//
//    auto results = eltC(C);
//    std::cout << "Result = " << results << std::endl;
