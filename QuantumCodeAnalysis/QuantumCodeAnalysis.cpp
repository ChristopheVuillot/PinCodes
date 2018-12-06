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
    typedef Givaro::Modular<double> Field;
    Field F(p); 
    Field::Element zero = F.zero;
    Field::Element one = F.one;

    // Reading the matrices from a file
    Field::Element * CCX, * CCZ, * MULTXZ;
    size_t mx, nx, mz, nz;
    // Specify the sparse format
    FFLAS::FFLAS_FORMAT format;
    format = FFLAS::FflasSMS;

    FFLAS::ReadMatrix(filex.c_str(), F, mx, nx, CCX, format);
    FFLAS::ReadMatrix(filez.c_str(), F, mz, nz, CCZ, format);
    MULTXZ = FFLAS::fflas_new(F,mx,mz);

    // Multiplying them
    FFLAS::FFLAS_TRANSPOSE tx, tz;
    tx = FFLAS::FflasNoTrans;
    tz = FFLAS::FflasTrans;
    FFLAS::fgemm(F, tx, tz, mx, mz, nx, one, CCX, nx, CCZ, nz, zero, MULTXZ, nz);

    // Result should be all zeros
    FFLAS::WriteMatrix (std::cout << "MULTXZ:=", F, mx, mz, MULTXZ, mz) << std::endl;
    // FFLAS::fflas_delete(CCX);
    // FFLAS::fflas_delete(CCZ);
    // FFLAS::fflas_delete(MULTXZ);

    return 0;
}
