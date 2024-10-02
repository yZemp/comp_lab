#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    printf("(.7 + .1) + .3\n");
    printf("%.16f\n", ((.7 + .1) + .3));

    printf("\n");

    printf(".7 + (.1 + .3)\n");
    printf("%.16f\n", (.7 + (.1 + .3)));

    printf("\n");


    float xt = 1.e20;
    float yt = -1.e20;
    float zt = 1;

    printf("Floats:\n");
    printf("xt: %f\tyt: %f\tzt: %f\n", xt, yt, zt);

    printf("(xt + yt) + zt\n");
    printf("%.16f\n", ((xt + yt) + zt));

    printf("\n");

    printf("xt + (yt + zt)\n");
    printf("%.16f\n", (xt + (yt + zt)));

    printf("\n");
    
    
    double xt1 = 1.e20;
    double yt1 = -1.e20;
    double zt1 = 1;

    printf("Doubles:\n");
    printf("xt: %f\tyt: %f\tzt: %f\n", xt1, yt1, zt1);

    printf("(xt + yt) + zt\n");
    printf("%.16f\n", ((xt1 + yt1) + zt1));

    printf("\n");

    printf("xt + (yt + zt)\n");
    printf("%.16f\n", (xt1 + (yt1 + zt1)));

    printf("\n");

    
    return(0);
}