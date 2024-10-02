#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    printf("First cycle\n");
    
    printf("Single:\n");
    float f = 1.2e34;
    for (size_t i = 0; i < 24; i++) {
        f *= 2;
        printf("%d:\t%f\n", i, f);
    }
    printf("Double:\n");
    double d = 1.2e34;
    for (size_t i = 0; i < 24; i++) {
        d *= 2;
        printf("%d:\t%f\n", i, d);
    }
    
    printf("Second cycle\n");

    printf("Signle:\n");
    float f = 1.2e304;
    for (size_t i = 0; i < 24; i++) {
        f /= 2;
        printf("%d:\t%f\t%f\n", i, f, 1 + f);
    }

    printf("Double:\n");
    double d = 1.2e304;
    for (size_t i = 0; i < 24; i++) {
        d /= 2;
        printf("%d:\t%f\t%f\n", i, d, 1 + d);
    }

    return(0);
}