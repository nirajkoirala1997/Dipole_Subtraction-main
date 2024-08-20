#include "ampvv.h" // Include the header file declaring your C++ function

extern "C" {
    void cpp_wrapper_(double *sn, double *un, double *m0m0, double *m0m1, double *m1m1, double *m0m2a, double *m0m2b, double *m0m2c, double *m0m2d) {
        // Call your C++ function from here
        // Assume your_cpp_function takes two double inputs and returns seven doubles

        ampvv(*sn, *un, *m0m0, *m0m1, *m1m1, *m0m2a, *m0m2b, *m0m2c, *m0m2d);

    }
}
