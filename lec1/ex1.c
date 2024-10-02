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
    fptr = fopen("data.txt", "w");


    for (size_t i = 0.1; i <= 1; i += .1) {

        fprintf(fptr, "Hello World");
    }
    
    // Close the file
    fclose(fptr); 

    return(0);
}