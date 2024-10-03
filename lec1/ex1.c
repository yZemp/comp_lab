#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define MAX_ITER 6 /* NUMERO MAX ITER SERIE*/

double exp_series(float x, int n);
double int_power(double x, int n);
double err_scale(double x, int n);
int fact(int n);

int main(int argc, char *argv[]) {
    int N = 20; // NUMERO CAMPIONAMENTI X
    double arr[N];
    double start = 0.1;
    double end = 1.0;
    double step = (end - start) / (N - 1); // Calculate step size for linear spacing

    // Fill an array with linearly spaced values from 0.1 to 1
    for (int i = 0; i < N; i++) {
        arr[i] = start + i * step;
    }

    double **data = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        data[i] = (double *) malloc(20 * sizeof(double));
    }

    // Calculate the matrix elements
    for (int i = 1; i < MAX_ITER; i++) {
        printf("Calculating series data: i = %d\n", i);
        for (int j = 0; j < N; j++) {
            data[i][j] = exp_series(arr[j], i);
        }
    }

    // Create a file
    FILE *fptr;
    fptr = fopen("data.dat", "w");
        
    FILE *ferr;
    ferr = fopen("error_scale.dat", "w");

    // Printo tutto l'array (prima riga)
    for (int i = 0; i < N; i++) {
        fprintf(fptr, "%lf ", arr[i]);
    }

    fprintf(fptr, "\n");

    // Printo tutto l'esponenziale (seconda riga)
    for (int i = 0; i < N; i++) {
        fprintf(fptr, "%lf ", exp(arr[i]));
    }

    fprintf(fptr, "\n");


    for (int i = 1; i < MAX_ITER; i++) {

        // Scrivo i valori della riga corrispondente della matrice
        for (int j = 0; j < N; j++) {
            fprintf(fptr, "%lf ", data[i][j]);
        }

        // Salvo anche l'errore che mi aspetto
        for (int j = 0; j < N; j++) {
            fprintf(ferr, "%lf ", err_scale(arr[j], i));
        }

        // A capo dopo ogni riga
        fprintf(fptr, "\n");
        fprintf(ferr, "\n");
    }

    // Close the file
    fclose(fptr); 
    fclose(ferr);

    return(0);
}


double exp_series(float x, int n) {
    double g = 0;
    for (size_t i = 0; i < n; i++) {
        g += int_power(x, i) / fact(i);
        // printf("%f\n", g);
    }
    
    return g;
}

double int_power(double x, int n) {
    if(n == 0) return 1;
    if(n == 1) return x;
    return x * int_power(x, n - 1);
}

int fact(int n) {
    if(n == 1 || n == 0) return 1;
    return n * fact(n - 1);
}

double err_scale(double x, int n) {
    return int_power(x, n + 1) / fact(n + 1);
}