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
    bool xperiodic = input.getInt("IsPeriodicX");
    bool yperiodic = input.getInt("IsPeriodicY");

    double Kx = input.getInt("Kx");
    double Ky = input.getInt("Ky");
    double Kz = input.getInt("Kz");

    double Hx = input.getInt("Hx");
    double Hy = input.getInt("Hy");
    double Hz = input.getInt("Hz");

    // sweep parameters
    auto totalSweeps = input.getInt("totalSweeps");
    auto sw_table = InputGroup(input,"sweep_table");
    auto sweeps = Sweeps(totalSweeps,sw_table);
    println(sweeps);

    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    auto ampo = AutoMPO(sites);
    std::cout << " Reading from File " << std::endl;
    readFromFile("sites.dat",sites);
    auto psi = readFromFile<MPS>(std::string("psi.dat"),sites);


    // =================== Observables =========================
    // ---------------------- On-site moments ------------------
    std::cout << " Starting calculations of Observables " << std::endl;
    std::ofstream outfile;
    outfile.open("Observables_onsite.dat");
    println("\nj Sx Sy Sz = ");

    for (int i=1; i<=N; i++) {
        psi.position(i);
        ITensor ket = psi.A(i);
        ITensor bra = dag(prime(ket,"Site"));

        ITensor Lxjop = sites.op("Sx",i); //*sites.op("Sx",j);
        ITensor Lyjop = sites.op("Sy",i); //*sites.op("Sy",j);
        ITensor Lzjop = sites.op("Sz",i); //*sites.op("Sz",j);

        //take an inner product
        auto Sxj = (bra*Lxjop*ket).cplx();
        auto Syj = (bra*Lyjop*ket).cplx();
        auto Szj = (bra*Lzjop*ket).cplx();

        printfln("%d %.12f %.12f %.12f",i,Sxj,Syj,Szj);
        outfile << i << " "  << Sxj << "  " << Syj << "  " << Szj << "  \n";
    }
    outfile.close();


    // ------------------- Correlation functions ---------------
    // Sx Sx correlators
    std::cout << "Starting calculations of Sx-Sx correlators:" << std::endl;
    std::ofstream outfileX;
    std::string opA = "Sx";
    std::string opB = "Sx";
    std::string outfilename="Observables_"+opA+opB+".dat";
    outfileX.open(outfilename);

    Eigen::MatrixXcd Obs(N, N); Obs.setZero();
    outfileX << std::setprecision(12);
    std::cout << std::setprecision(12);

    outfileX << opA << opB << " = (" << N << "," << N << ") \n";
    std::cout << opA << opB << " = (" << N << "," << N << ") \n";

    for (int i=1; i<=N; i++) {
        psi.position(i);
        ITensor op_i = sites.op(opA, i);

        for (int j=1; j<=N; j++) {
            if (j<=i) {
                outfileX << "(0,0)" << " ";
                std::cout << "(0,0)" << " ";
            } else {
                ITensor op_j = op(sites, opA, j);

                ITensor C = psi(i);
                C *= op_i;
                auto ir = commonIndex(psi(i),psi(i+1),"Link");
                C *= dag(prime(prime(psi(i),"Site"),ir));

                for(int k = i+1; k < j; ++k) {
                    C *= psi(k);
                    C *= dag(prime(psi.A(k),"Link"));
                }

                C *= psi(j);
                C *= op_j;

                auto jl = commonIndex(psi(j),psi(j-1),"Link"); //index linking j to j-1:
                C *= dag(prime(prime(psi(j),"Site"), jl));


                auto result = eltC(C);
                outfileX << result << " ";
                std::cout << result << " ";
            }
        }
        outfileX << " \n";
        std::cout << " \n";
    }
    outfileX.close();

    // Sy Sy correlators
    std::cout << "Starting calculations of Sx-Sx correlators:" << std::endl;
    std::ofstream outfileY;
    opA = "Sy";
    opB = "Sy";
    outfilename="Observables_"+opA+opB+".dat";
    outfileY.open(outfilename);

    Obs.setZero();
    outfileY << std::setprecision(12);
    std::cout << std::setprecision(12);

    outfileY << opA << opB << " = (" << N << "," << N << ") \n";
    std::cout << opA << opB << " = (" << N << "," << N << ") \n";

    for (int i=1; i<=N; i++) {
        psi.position(i);
        ITensor op_i = sites.op(opA, i);

        for (int j=1; j<=N; j++) {
            if (j<=i) {
                outfileY << "(0,0)" << " ";
                std::cout << "(0,0)" << " ";
            } else {
                ITensor op_j = op(sites, opA, j);

                ITensor C = psi(i);
                C *= op_i;
                auto ir = commonIndex(psi(i),psi(i+1),"Link");
                C *= dag(prime(prime(psi(i),"Site"),ir));

                for(int k = i+1; k < j; ++k) {
                    C *= psi(k);
                    C *= dag(prime(psi.A(k),"Link"));
                }

                C *= psi(j);
                C *= op_j;

                auto jl = commonIndex(psi(j),psi(j-1),"Link"); //index linking j to j-1:
                C *= dag(prime(prime(psi(j),"Site"), jl));


                auto result = eltC(C);
                outfileY << result << " ";
                std::cout << result << " ";
            }
        }
        outfileY << " \n";
        std::cout << " \n";
    }
    outfileY.close();

    // Sz Sz correlators
    std::cout << "Starting calculations of Sx-Sx correlators:" << std::endl;
    std::ofstream outfileZ;
    opA = "Sz";
    opB = "Sz";
    outfilename="Observables_"+opA+opB+".dat";
    outfileZ.open(outfilename);

    Obs.setZero();
    outfileZ << std::setprecision(12);
    std::cout << std::setprecision(12);

    outfileZ << opA << opB << " = (" << N << "," << N << ") \n";
    std::cout << opA << opB << " = (" << N << "," << N << ") \n";

    for (int i=1; i<=N; i++) {
        psi.position(i);
        ITensor op_i = sites.op(opA, i);

        for (int j=1; j<=N; j++) {
            if (j<=i) {
                outfileZ << "(0,0)" << " ";
                std::cout << "(0,0)" << " ";
            } else {
                ITensor op_j = op(sites, opA, j);

                ITensor C = psi(i);
                C *= op_i;
                auto ir = commonIndex(psi(i),psi(i+1),"Link");
                C *= dag(prime(prime(psi(i),"Site"),ir));

                for(int k = i+1; k < j; ++k) {
                    C *= psi(k);
                    C *= dag(prime(psi.A(k),"Link"));
                }

                C *= psi(j);
                C *= op_j;

                auto jl = commonIndex(psi(j),psi(j-1),"Link"); //index linking j to j-1:
                C *= dag(prime(prime(psi(j),"Site"), jl));


                auto result = eltC(C);
                outfileZ << result << " ";
                std::cout << result << " ";
            }
        }
        outfileZ << " \n";
        std::cout << " \n";
    }
    outfileY.close();

    return 0;
}
