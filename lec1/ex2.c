#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    printf("(.7 + .1) + .3\n");
    printf("%.16e\n", ((.7 + .1) + .3));

    printf("\n");

    printf(".7 + (.1 + .3)\n");
    printf("%.16e\n", (.7 + (.1 + .3)));

    printf("\n");

    float xt = 1.e20;
    float yt = -1.e20;
    float zt = 1;

    printf("(xt + yt) + zt\n");
    printf("%.16e\n", ((xt + yt) + zt));

    printf("\n");

    printf("xt + (yt + zt)\n");
    printf("%.16e\n", (xt + (yt + zt)));

    printf("\n");

    return(0);
}