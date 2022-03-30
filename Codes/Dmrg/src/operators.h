//
// Created by shifeng on 3/30/22.
//

#ifndef DMRG_OPERATORS_H
#define DMRG_OPERATORS_H
typedef std::complex<double> dcomplex;
namespace itensor {
    ITensor myOp (Index s, std::string name) {
        auto sp = prime(s);
        auto res = ITensor(s, sp);

        if (name == "Szp") {
            res.set(s=1, sp=1, +1.0);
        }
        else if (name == "Szm") {
            res.set(s=2, sp=2, +1.0);
        }
        else if (name == "Sxp") {
            res.set(s=1, sp=1, +0.5);
            res.set(s=1, sp=2, +0.5);
            res.set(s=2, sp=1, +0.5);
            res.set(s=2, sp=2, +0.5);
        }
        else if (name == "Sxm") {
            res.set(s=1, sp=1, +0.5);
            res.set(s=1, sp=2, -0.5);
            res.set(s=2, sp=1, -0.5);
            res.set(s=2, sp=2, +0.5);
        }
        else if (name == "Syp") {
            res.set(s=1, sp=1, +0.5);
            res.set(s=1, sp=2, dcomplex(0, 0.5));
            res.set(s=2, sp=1, dcomplex(0, 0.5));
            res.set(s=2, sp=2, -0.5);
        }
        else if (name == "Sym") {
            res.set(s=1, sp=1, +0.5);
            res.set(s=1, sp=2, dcomplex(0, -0.5));
            res.set(s=2, sp=1, dcomplex(0, -0.5));
            res.set(s=2, sp=2, -0.5);
        }
        else {
            std::cout << name << std::endl;
            Error("Operator name not defined");
        }

        return res;
    }
}
#endif //DMRG_OPERATORS_H
