#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    printf("First cycle\n");

    printf("Single:\n");
    float f1 = 1.2e34;
    for (size_t i = 0; i < 24; i++) {
        f1 *= 2;
        printf("%d:\t%.20e\n", i, f1);
    }
    printf("Double:\n");
    double d1 = 1.2e304;
    for (size_t i = 0; i < 24; i++) {
        d1 *= 2;
        printf("%d:\t%.20e\n", i, d1);
    }
    
    printf("Second cycle\n");

    printf("Single:\n");
    float f2 = 1e-13;
    for (size_t i = 0; i < 24; i++) {
        f2 /= 2;
        printf("%d:\t%.20e\t%.20e\n", i, f2, 1 + f2);
    }

    printf("Double:\n");
    double d2 = 1e-13;
    for (size_t i = 0; i < 24; i++) {
        d2 /= 2;
        printf("%d:\t%.20e\t%.20e\n", i, d2, 1 + d2);
    }

    return(0);
}