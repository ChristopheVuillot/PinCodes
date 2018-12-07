#include <givaro/modular.h>
#include <iostream>

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/fflas_io.h"

int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "Usage: analyse <p> <ccx> <ccz>" << std::endl;
        return -1;
    }

    int p = atoi(argv[1]);
    std::string filex = argv[2];
    std::string filez = argv[3];

    // Creating the finite field Z/pZ
    typedef Givaro::Modular<float> Field;
    Field F(p); 
    Field::Element zero = F.zero;
    Field::Element one = F.one;

    // Reading the matrices from a file
    Field::Element * CCX, * CCZ, * MULTXZ;
    size_t mx, nx, mz, nz;
    // Specify the sparse format
    FFLAS::FFLAS_FORMAT format;
    format = FFLAS::FflasSMS;

    std::cout << "Reading matrices from file: " << std::endl;
    FFLAS::ReadMatrix(filex.c_str(), F, mx, nx, CCX, format);
    FFLAS::ReadMatrix(filez.c_str(), F, mz, nz, CCZ, format);
    MULTXZ = FFLAS::fflas_new(F,mx,mz);

    // Multiplying them
    FFLAS::FFLAS_TRANSPOSE tx, tz;
    tx = FFLAS::FflasNoTrans;
    tz = FFLAS::FflasTrans;
    std::cout << "Multiplying matrices: " << std::endl;
    FFLAS::fgemm(F, tx, tz, mx, mz, nx, one, CCX, nx, CCZ, nz, zero, MULTXZ, nz);

    // Result should be all zeros
    Field::Element * onesL, * onesR, * rowsum;
    onesL = FFLAS::fflas_new(F,mx,1);
    rowsum = FFLAS::fflas_new(F,mx,1);
    onesR = FFLAS::fflas_new(F,mz,1);
    for (size_t i=0; i<mx; i++) {
	    F.init(onesL[i], one);
	    F.init(rowsum[i], zero);
    }
    for (size_t i=0; i<mz; i++) {
	    F.init(onesR[i], one);
    }


    std::cout << "Checking if commute: " << std::endl;
    FFLAS::fgemv(F, FFLAS::FflasNoTrans, mx, mz, one, MULTXZ, mz, onesR, 1, zero, rowsum, 1);
    std::cout << "Checking if commute: " << std::endl;
    Field::Element res = FFLAS::fdot(F, mx, onesL, 1, rowsum, 1);

    std::cout << "HxHz^T == 0 ? -> " << res << std::endl;

    // FFLAS::WriteMatrix (std::cout << "MULTXZ:=", F, mx, mz, MULTXZ, mz) << std::endl;
    std::cout << "plop" << std::endl;
    FFLAS::fflas_delete(CCX);
    std::cout << "plop" << std::endl;
    // FFLAS::fflas_delete(CCZ);
    std::cout << "plop" << std::endl;
    FFLAS::fflas_delete(MULTXZ);
    std::cout << "plop" << std::endl;
    FFLAS::fflas_delete(onesL);
    std::cout << "plop" << std::endl;
    FFLAS::fflas_delete(onesR);
    std::cout << "plop" << std::endl;
    FFLAS::fflas_delete(rowsum);
    std::cout << "plop" << std::endl;

    return 0;
}
