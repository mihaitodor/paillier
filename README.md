# paillier

A self-contained implementation of the Paillier cryptosystem with the CRT (Chinese Remainder Theorem) optimisation for decryption using GMP.

[SeComLib](https://github.com/mihaitodor/SeComLib) was used as a reference and this implementation contains the same optimisations. For further details, please consult the following [code](https://github.com/mihaitodor/SeComLib/blob/master/core/paillier.cpp) and paper: "Public-Key Cryptosystems Based on Composite Degree Residuosity Classes" by Pacal Paillier, 1999.

To build and run the code: `gcc -lgmp paillier.cpp && ./a.out 1 2`
