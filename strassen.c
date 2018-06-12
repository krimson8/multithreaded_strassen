#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>


int n;

int NUM_THREAD;

typedef struct _Matrix {
    int **v;
    int size;
} Matrix;

int thread_count = 4;
pthread_mutex_t mutexA = PTHREAD_MUTEX_INITIALIZER;
int counter = 0;

void matrix_create(Matrix *thiz, int size, FILE *fp)
{
    assert(thiz != NULL);   
    thiz->size = size;

    // dynamic allocate 2d matrix
    int *arr = (int*)malloc(size * size * sizeof(int));
    thiz->v = (int**)malloc(size * sizeof(int*));
    for (int i = 0; i < size; ++i)
        thiz->v[i] = arr + i*size;

    if (fp) {
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                fscanf(fp, "%d", &(thiz->v[i][j]));
    }
}

void matrix_read(Matrix *thiz, FILE *fp)
{
    int m, n;
    fscanf(fp, "%d %d", &m, &n);
    assert(m == n);
    matrix_create(thiz, m, fp);
}

void matrix_print(Matrix thiz)
{
    for (int i = 0; i < thiz.size; ++i) {
        for (int j = 0; j < thiz.size; ++j)
            printf("%d ", thiz.v[i][j]);
        puts("");
    }
    puts("");
}

void matrix_try_create(Matrix *thiz, int size, FILE *fp)
{
    if (thiz->size == 0)
        matrix_create(thiz, size, NULL);
}

void matrix_check_AB_try_create_C(const Matrix *A, const Matrix *B, Matrix *C)
{
    assert(A != NULL && B != NULL);
    assert(A->size == B->size);
    matrix_try_create(C, A->size, NULL);
}

void matrix_add(const Matrix *A, const Matrix *B, Matrix *C)
{
    matrix_check_AB_try_create_C(A, B, C);
    int size = A->size;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            C->v[i][j] = A->v[i][j] + B->v[i][j];
        }
    } 
}

void matrix_sub(const Matrix *A, const Matrix *B, Matrix *C)
{
    matrix_check_AB_try_create_C(A, B, C);
    int size = A->size;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            C->v[i][j] = A->v[i][j] - B->v[i][j];
        }
    } 
}

void matrix_mul(Matrix *A, Matrix *B, Matrix *C)
{
    matrix_check_AB_try_create_C(A, B, C);
    int size = A->size;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            C->v[i][j] = 0;
            for (int k = 0; k < size; ++k) {
                C->v[i][j] += A->v[i][k] * B->v[k][j];
            }
        }
    } 
}

void matrix_divide_4(Matrix *thiz, Matrix *block)
{
    int size = thiz->size;
    int size_divided = size/2;
    for (int i = 0; i < 4; ++i)
        matrix_create(&block[i], size_divided, 0);
    for (int j = 0; j < size_divided; ++j) {
        block[0].v[j] = thiz->v[j]; // TODO double allocate
        block[1].v[j] = thiz->v[j] + size_divided;
        block[2].v[j] = thiz->v[j+size_divided];
        block[3].v[j] = thiz->v[j+size_divided] + size_divided;
    }
}

void matrix_combine_4(Matrix *thiz, Matrix *block)
{
    int size_divided = block[0].size;
    int size= 2*size_divided;
    matrix_try_create(thiz, size, NULL);
    for (int i = 0; i < size_divided; ++i) {
        for (int j = 0; j < size_divided; ++j) {
            thiz->v[i][j] = block[0].v[i][j];
            thiz->v[i][j+size_divided] = block[1].v[i][j];
            thiz->v[i+size_divided][j] = block[2].v[i][j];
            thiz->v[i+size_divided][j+size_divided] = block[3].v[i][j];
        }
    }
}

void _strassen_mul(const Matrix *A, const Matrix *B, Matrix **M, Matrix **tmp)

void strassen_mul(Matrix *A_all, Matrix *B_all, Matrix *C_all, bool parallel, int id)
{
    Matrix A[4] = {0}, 
           B[4] = {0}, 
           C[4] = {0}, 
           M[7] = {0}, 
           tmp[10] = {0};

    matrix_divide_4(A_all, A);
    matrix_divide_4(B_all, B);

    _strassen_mul(A, B, &M, &tmp, id);

    //********** parallel if needed **************
    if (!parallel || id == 0) {
        matrix_add(&A[0], &A[3], &tmp[0]);
        matrix_add(&B[0], &B[3], &tmp[1]);
        matrix_mul(&tmp[0], &tmp[1], &M[0]);
    }

    if (!parallel || id == 1) {
        matrix_add(&A[2], &A[3], &tmp[2]);
        matrix_mul(&tmp[2], &B[0], &M[1]);
    }

    if (!parallel || id == 2) {
        matrix_sub(&B[1], &B[3], &tmp[3]);
        matrix_mul(&A[0], &tmp[3], &M[2]);
    }

    if (!parallel || id == 3) {
        matrix_sub(&B[2], &B[0], &tmp[4]);
        matrix_mul(&A[3], &tmp[4], &M[3]);
    }


    if (!parallel || id == 4) {
        matrix_add(&A[0], &A[1], &tmp[5]);
        matrix_mul(&tmp[5], &B[3], &M[4]);
    }

    if (!parallel || id == 5) {
        matrix_sub(&A[2], &A[0], &tmp[6]);
        matrix_add(&B[0], &B[1], &tmp[7]);
        matrix_mul(&tmp[6], &tmp[7], &M[5]);
    }

    if (!parallel || id == 6) {
        matrix_sub(&A[1], &A[3], &tmp[8]);
        matrix_add(&B[2], &B[3], &tmp[9]);
        matrix_mul(&tmp[8], &tmp[9], &M[6]);
    }

    //********** end of parallel section **************

    matrix_add(&M[0], &M[3], &tmp[0]);
    matrix_sub(&tmp[0], &M[4], &tmp[0]);
    matrix_add(&tmp[0], &M[6], &C[0]);
    matrix_add(&M[2], &M[4], &C[1]);
    matrix_add(&M[1], &M[3], &C[2]);
    matrix_sub(&M[0], &M[1], &tmp[0]);
    matrix_add(&tmp[0], &M[2], &tmp[0]);
    matrix_add(&tmp[0], &M[5], &C[3]);

    matrix_combine_4(C_all, C);
}

void* child(void *arg)
{
    pthread_exit(NULL);
}

int main(int argc, const char **argv) 
{ 

    if (argc != 4) {
        puts("./main path/to/test/data type_of_mul num_thread");
        exit(EXIT_FAILURE);
    }

    FILE *fp = fopen(argv[1], "r");
    if (fp==NULL) {
        perror(argv[1]);
        exit(EXIT_FAILURE);
    }

    NUM_THREAD = atoi(argv[3]);

    Matrix A = {0}, 
           B = {0}, 
           C = {0};
    memset(&A, 0, sizeof(Matrix));
    memset(&B, 0, sizeof(Matrix));
    memset(&C, 0, sizeof(Matrix));
    matrix_read(&A, fp);
    matrix_read(&B, fp);
    fclose(fp);

    clock_t begin, end; 
    double time_spent;

    begin = clock();
    if (atoi(argv[2]) == 0) {
        matrix_mul(&A, &B, &C);
    } else {
        strassen(&A, &B, &C);
    }
    end = clock();
    time_spent = (double)(end - begin) / (CLOCKS_PER_SEC);
    printf("%lf\n",time_spent);    

    return 0;



    pthread_t t[thread_count];
    int param_i[thread_count];
    for(int i = 0; i < thread_count; i++) {
        param_i[i] = i;
        pthread_create(&t[i], NULL, child, &param_i[i]);
    }

    for(int i = 0; i < thread_count; i++) {
        pthread_join(t[i], NULL);
    }


}
