# This is testbot
# The code in this file runs anytime code is pushed or a pull request is made to the main branch

println!("test")

using Oscar
using Pkg
Pkg.add(OscarPuiseuxPolynomial)
using OscarPuiseuxPolynomial

Kt, (t1,t2) = puiseux_polynomial_ring(QQ, ["t1","t2"])
