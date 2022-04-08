#include <iostream>
#include "itensor/all.h"
#include "operators.h"
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
    int spin2pauli = input.getInt("spin2pauli");

    // ================= Read in mask ====================================
    int maskSize = input.getInt("maskSize");
    std::vector<int> mpsIndx(maskSize);
    std::vector<int> x(maskSize);
    std::vector<int> y(maskSize);
    std::vector<int> z(maskSize);
    std::vector<std::string> maskList(maskSize);
    auto mask = InputGroup(input, "mask");
    printfln("mpsId x y z  maskname");
    if (mask.GotoGroup()) {
        mask.SkipLine();
        for (int i = 0; i<maskSize; i++) {
            mask.file() >> mpsIndx[i] >> x[i] >> y[i] >> z[i];
            assert(x[i] + y[i] + z[i] == 1);
            if (x[i] == 1) {maskList[i] = "Sx+";}
            else if (x[i] == -1) {maskList[i] = "Sx-";}
            else {
                if (y[i] == 1) {maskList[i] = "Sy+";}
                else if (y[i] == -1) {maskList[i] = "Sy-";}
                else {
                    if (z[i] == 1) {maskList[i] = "Sz+";}
                    else if (z[i] == -1) {maskList[i] = "Sz-";}
                }
            }
            printfln("  %d   %d %d %d   %d",mpsIndx[i],x[i],y[i],z[i],maskList[i]);
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


    // ============== Build RDM from MPS =================
    int start = *std::min_element(mpsIndx.begin(), mpsIndx.end()); println("Mask start at site: ", start);
    int end = *std::max_element(mpsIndx.begin(), mpsIndx.end()); println("Mask ends at site: ", end);

    psi.position(start);
    auto psidag = dag(psi);
    psidag.prime("Link");

    auto li_1 = leftLinkIndex(psi, start);
    auto rho = prime(psi(start), li_1)*prime(psidag(start), "Site");
    // the bulk of tensor train
    for (int i = start+1; i < end; i++) { // i-th neighbor to the right of the head
        if (std::find(mpsIndx.begin(), mpsIndx.end(), i) != mpsIndx.end()) {
            // if i+1 is in mpsIndx (continuous indices)
            println("if true");
            auto li = rightLinkIndex(psi, i);
            rho *= prime(psi(i), li);
            rho *= prime(psidag(i), "Site");
        } else {
            println("else");
            rho *= psi(i);
            rho *= psidag(i);
        }
    }

//    auto lend = rightLinkIndex(psi, end);
//    rho *= prime(psi(end), lend);
//    rho *= prime(psidag(end), "Site");





//    // =================== Multi-Point Observables =========================
//    int start = *std::min_element(mpsIndx.begin(), mpsIndx.end()); println("Mask start at site: ", start);
//    int end = *std::max_element(mpsIndx.begin(), mpsIndx.end()); println("Mask ends at site: ", end);
//
//    // initialize the tensor train from the left side
//    psi.position(start);
//    ITensor C = psi(start);
//
//    // the head of tensor train
//    std::string S_op;
//    S_op = maskList[0];
//    println("\nhead op = ", S_op, " at site ", start);
//    C *= myOp(sites(start), S_op); // op(sites, S_op, start);
//    auto ir = commonIndex(psi(start),psi(start+1),"Link");
//    C *= dag(prime(prime(psi(start), "Site"), ir));
//    C *= psi(start+1);
//
//    // the bulk of tensor train
//    int countId = 0;
//    for (int i = start+1; i < end; i++) { // i-th neighbor to the right of the head
//        if (std::find(mpsIndx.begin(), mpsIndx.end(), i) != mpsIndx.end()) {
//            // if i+1 is in mpsIndx (continuous indices)
//            S_op = maskList[i-start-countId];
//            println("bulk op = ", S_op, " at site ", i);
//            C *= myOp(sites(i), S_op); //op(sites, S_op, i);
//            C *= dag(prime(prime(psi(i), "Site"), "Link"));
//            C *= psi(i+1);
//        } else {
//            println("directly contract trivial bulk tensor at ", i);
//            countId ++;
//            C *= dag(prime(psi(i), "Link"));
//            C *= psi(i+1);
//        }
//    }
//
//    // the tail of tensor train
//    S_op = maskList[maskSize-1];
//    println("tail op = ", S_op, " at site ", end, '\n');
//    C *= myOp(sites(end), S_op);
//    auto il = commonIndex(psi(end-1),psi(end),"Link");
//    C *= dag(prime(prime(psi(end), "Site"), il));
//
//    auto results = eltC(C);
//    if (spin2pauli)
//        results *= pow(2,maskSize);
//    std::cout << "Result = " << results << std::endl;

//        // ====test====
//        ITensor Sz_2 = op(sites, "Sx", 2);
//        ITensor Sz_3 = op(sites, "Sy", 3);
//        ITensor Sz_4 = op(sites, "Sz", 4);
//        ITensor Sz_7 = op(sites, "Sz", 7);
//        ITensor Sz_8 = op(sites, "Sy", 8);
//        ITensor Sz_9 = op(sites, "Sx", 9);
//
//        psi.position(2);
//
//        ITensor C = psi(2);
//        C *= Sz_2;  // with primed physical index
//        auto ir = commonIndex(psi(2),psi(3),"Link");
//        C *= dag(prime(prime(psi(2), "Site"), ir));
//
//        C *= psi(3);
//        C *= Sz_3;
//        C *= dag(prime(prime(psi(3), "Site"), "Link"));
//
//        C *= psi(4);
//        C *= Sz_4;
//        C *= dag(prime(prime(psi(4), "Site"), "Link"));
//
//        C *= psi(5);
//        C *= dag(prime(psi(5), "Link"));
//
//        C *= psi(6);
//        C *= dag(prime(psi(6), "Link"));
//
//        C *= psi(7);
//        C *= Sz_7;
//        C *= dag(prime(prime(psi(7), "Site"), "Link"));
//
//        C *= psi(8);
//        C *= Sz_8;
//        C *= dag(prime(prime(psi(8), "Site"), "Link"));
//
//        C *= psi(9);
//        C *= Sz_9;
//        auto il = commonIndex(psi(8),psi(9),"Link");
//        C *= dag(prime(prime(psi(9), "Site"), il));
//
//        auto results = eltC(C * 64);
//        std::cout << "Result = " << results << std::endl;

    return 0;
}
