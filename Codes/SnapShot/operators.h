//
// Created by shifeng on 3/30/22.
//

#ifndef OBSERVE_OPERATORS_H
#define OBSERVE_OPERATORS_H
typedef std::complex<double> dcomplex;
namespace itensor {
    ITensor myOp (Index s, std::string name) {
        auto sp = prime(s);
        auto res = ITensor(s, sp);

        if (name == "Sz+") {
            res.set(s=1, sp=1, +1.0);
        }
        else if (name == "Sz-") {
            res.set(s=2, sp=2, +1.0);
        }
        else if (name == "Sx+") {
            res.set(s=1, sp=1, +0.5);
            res.set(s=1, sp=2, +0.5);
            res.set(s=2, sp=1, +0.5);
            res.set(s=2, sp=2, +0.5);
        }
        else if (name == "Sx-") {
            res.set(s=1, sp=1, +0.5);
            res.set(s=1, sp=2, -0.5);
            res.set(s=2, sp=1, -0.5);
            res.set(s=2, sp=2, +0.5);
        }
        else if (name == "Sy+") {
            res.set(s=1, sp=1, +0.5);
            res.set(s=1, sp=2, dcomplex(0, 0.5));
            res.set(s=2, sp=1, dcomplex(0, 0.5));
            res.set(s=2, sp=2, -0.5);
        }
        else if (name == "Sy-") {
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
#endif //OBSERVE_OPERATORS_H
