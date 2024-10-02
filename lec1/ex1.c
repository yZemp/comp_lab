#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double exp_series(float x) {
    // TODO
    return 0;
}

int main(int argc, char *argv[]) {
    // FILE *fptr;

    // Create a file
    // fptr = fopen("data.txt", "w");


    for (float x = 0.1; x <= 1; x += .1) {
        printf("%f\t%f\n", x, exp(x));
        // fprintf(fptr, "Hello World");
    }
    
    // Close the file
    // fclose(fptr); 

    return(0);
}