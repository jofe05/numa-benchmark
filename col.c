#include <stdio.h>
#include <stdlib.h>

int main(){

    size_t i,j;

    size_t dim = 81600 ;
    typedef int **Matrix;

    Matrix m = malloc(sizeof(int *) * dim);

	for ( i = 0; i < dim; i++ )
    		m[i] = malloc(sizeof(int) * dim);



    for (i = 0;i< dim; i++){
        for (j=0;j <dim;j++){
            m[j][i]= i*j;
        }
    }
    return 0;
}
