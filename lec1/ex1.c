#include <stdlib.h>
#include <stdio.h>
#include <math.h>


double exp_series(float x, int n);
double int_power(double x, int n);
int fact(int n);

int main(int argc, char *argv[]) {
    FILE *fptr1;
    FILE *fptr2;

    // Create a file
    fptr1 = fopen("exp.dat", "w");
    fptr2 = fopen("series.dat", "w");


    for (float x = 0.1; x <= 1; x += .05) {
        char* string1;
        if(0 > asprintf(&string1, "%f\t%f\n", x, exp(x))) return -1;
        fprintf(fptr1, string1);
        free(string1);

        char* string2;
        if(0 > asprintf(&string2, "%f\t%f\n", x, exp_series(x, 3))) return -1;
        fprintf(fptr2, string2);
        free(string2);

        // printf("%f\t%f\n", x, exp(x));
    }
    
    // Close the file
    fclose(fptr1); 
    fclose(fptr2); 

    return(0);
}


double exp_series(float x, int n) {
    double g = 0;
    for (size_t i = 1; i < n; i++) {
        g += int_power(x, i) / fact(i);
    }
    
    return g;
}

double int_power(double x, int n) {
    if(n == 1) return x;
    return x * int_power(x, n - 1);
}

int fact(int n) {
    if(n == 1) return n;
    return n * fact(n - 1);
}