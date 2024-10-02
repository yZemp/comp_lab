#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double exp_series(float x) {
    // TODO
    return 0;
}

int main(int argc, char *argv[]) {
    FILE *fptr;

    // Create a file
    fptr = fopen("data.dat", "w");


    for (float x = 0.1; x <= 1; x += .05) {
        char* string;
        if(0 > asprintf(&string, "%f\t%f\n", x, exp(x))) return -1;
        fprintf(fptr, string);
        free(string);

        printf("%f\t%f\n", x, exp(x));
    }
    
    // Close the file
    fclose(fptr); 

    return(0);
}