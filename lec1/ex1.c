#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    float f = 1.2e34;
    for (size_t i = 0; i < 24; i++) {
        f *= 2;
        printf("%d:\t%d", i, f);
    }
    

    return(0);
}