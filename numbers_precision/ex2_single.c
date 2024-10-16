#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define MAX 50e3 /* NUMERO MAX ITER PARTIAL SUM*/
#define PI 3.14159265358 /* Sono sicuro esista un metodo migliore di questo */

float series(int N);
float series_reversed(int N);
float int_power(float x, int n);

int main(int argc, char *argv[]) {

    // printf("Normal ordering series:\t\t%lf\n", series(MAX));
    // printf("Reversed ordering series:\t%lf\n", series_reversed(MAX));
    // printf("True value:\t\t\t%lf\n", int_power(PI, 2) / 6);

    // Create a file
    FILE *fnorm;
    fnorm = fopen("normal.dat", "w");
        
    FILE *frev;
    frev = fopen("reversed.dat", "w");

    float true_val = int_power(PI, 2) / 6;
    
    for (int i = 1; i < MAX; i += 1) {        
        if (!(i % 100)) {
            printf("Calculating partial sums: i = %d\n", i);
        }
        // Si potrebbe ottimizzare tenendo in memoria somma parziale a i - 1
        fprintf(fnorm, "%d %f\n", i, true_val - series(i));
        fprintf(frev, "%d %f\n", i, true_val - series_reversed(i));
    }
    
    // fprintf(fnorm, "\n");
    // fprintf(frev, "\n");

    // Close the file
    fclose(fnorm); 
    fclose(frev);

    return(0);
}


float series(int N) {
    float partial_sum = 0;
    for (int i = 1; i <= N; i++) {
        partial_sum += 1 / int_power(i, 2);
    }
    
    return partial_sum;
}

float series_reversed(int N) {
    float partial_sum = 0;
    for (int i = N; i > 0; i--) {
        partial_sum += 1 / int_power(i, 2);
    }
    
    return partial_sum;
}

float int_power(float x, int n) {
    if(n == 0) return 1;
    if(n == 1) return x;
    return x * int_power(x, n - 1);
}