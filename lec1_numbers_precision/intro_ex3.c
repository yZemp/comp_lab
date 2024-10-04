#include <stdlib.h>
#include <stdio.h>

// NOTE:
// Program frequently killed
// Not enough swap memory?

int main(int argc, char *argv[]) {
    float increment = .1;

    do {
        printf("Difference:\t%a\n", increment);
        printf("Difference:\t%.*f\n", increment);
        increment /= 2;
    } while (increment > 0);

    printf("\n");


    
    float epsilon = .1;
    float temp = .1;
    
    while ((1 + temp) != 1) {
        epsilon = temp;
        temp /= 2;
    }

    printf("Machine Epsilon (single):\t%a\n", epsilon);
    printf("Machine Epsilon (single):\t%.16f\n", epsilon);
    printf("Machine Epsilon (single):\t%.e\n", epsilon);
    printf("Single finished");
    
    

    double depsilon = .1;
    double dtemp = .1;
    
    while ((1 + dtemp) != 1) {
        depsilon = dtemp;
        dtemp /= 2;
    }

    printf("Machine Epsilon (double):\t%a\n", depsilon);
    printf("Machine Epsilon (double):\t%.16f\n", depsilon);
    printf("Machine Epsilon (double):\t%.e\n", depsilon);
    
    

    return(0);
}