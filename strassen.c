#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>


int n;

typedef struct _Matrix {
    int **v;
    int size;
} Matrix;

typedef struct _ThreadParam {
    Matrix *A, *B, *M, *tmp;
    int id;
} ThreadParam;

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

void matrix_mul(const Matrix *A, const Matrix *B, Matrix *C)
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

void _strassen_mul(const Matrix *A, const Matrix *B, Matrix *M, Matrix *tmp, const int id)
{
    if (id == -1 || id == 0) {
        matrix_add(&A[0], &A[3], &tmp[0]);
        matrix_add(&B[0], &B[3], &tmp[1]);
        matrix_mul(&tmp[0], &tmp[1], &M[0]);
    }

    if (id == -1 || id == 1) {
        matrix_add(&A[2], &A[3], &tmp[2]);
        matrix_mul(&tmp[2], &B[0], &M[1]);
    }

    if (id == -1 || id == 2) {
        matrix_sub(&B[1], &B[3], &tmp[3]);
        matrix_mul(&A[0], &tmp[3], &M[2]);
    }

    if (id == -1 || id == 3) {
        matrix_sub(&B[2], &B[0], &tmp[4]);
        matrix_mul(&A[3], &tmp[4], &M[3]);
    }


    if (id == -1 || id == 4) {
        matrix_add(&A[0], &A[1], &tmp[5]);
        matrix_mul(&tmp[5], &B[3], &M[4]);
    }

    if (id == -1 || id == 5) {
        matrix_sub(&A[2], &A[0], &tmp[6]);
        matrix_add(&B[0], &B[1], &tmp[7]);
        matrix_mul(&tmp[6], &tmp[7], &M[5]);
    }

    if (id == -1 || id == 6) {
        matrix_sub(&A[1], &A[3], &tmp[8]);
        matrix_add(&B[2], &B[3], &tmp[9]);
        matrix_mul(&tmp[8], &tmp[9], &M[6]);
    }
}

void* pthread_func(void *arg)
{
    ThreadParam *t = (ThreadParam*)arg; 
    _strassen_mul(t->A, t->B, t->M, t->tmp, t->id);
    pthread_exit(NULL);
}

void strassen_mul(Matrix *A_all, Matrix *B_all, Matrix *C_all, bool parallel)
{
    Matrix A[4] = {0}, 
           B[4] = {0}, 
           C[4] = {0}, 
           M[7] = {0}, 
           tmp[10] = {0};

    matrix_divide_4(A_all, A);
    matrix_divide_4(B_all, B);

    if (parallel) {
        int num_thread = 7;
        pthread_t thread[num_thread];
        ThreadParam thread_param[num_thread];

        for(int i = 0; i < num_thread; i++) {
            ThreadParam t = {.A=A, .B=B, .M=M, .tmp=tmp, .id = i};
            memcpy(&thread_param[i], &t, sizeof(ThreadParam));
            pthread_create(&thread[i], NULL, pthread_func, &thread_param[i]);
        }
        for(int i = 0; i < num_thread; i++) {
            pthread_join(thread[i], NULL);
        }
    } else {
        _strassen_mul(A, B, M, tmp, -1);
    }

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

int main(int argc, const char **argv) 
{ 

    if (argc != 3) {
        puts("./main path/to/test/data type_of_mul");
        exit(EXIT_FAILURE);
    }

    FILE *fp = fopen(argv[1], "r");
    if (fp==NULL) {
        perror(argv[1]);
        exit(EXIT_FAILURE);
    }
    int mul_type = atoi(argv[2]);

    Matrix A = {0}, 
           B = {0}, 
           C = {0};
    memset(&A, 0, sizeof(Matrix));
    memset(&B, 0, sizeof(Matrix));
    memset(&C, 0, sizeof(Matrix));
    matrix_read(&A, fp);
    matrix_read(&B, fp);
    fclose(fp);

    struct timeval start, end;
    gettimeofday(&start, NULL);

    switch (mul_type) {
    case 0:
        matrix_mul(&A, &B, &C);
        break;
    case 1:
        strassen_mul(&A, &B, &C, false);
        break;
    case 2:
        strassen_mul(&A, &B, &C, true);
        break;
    }

    gettimeofday(&end, NULL);
    int time_spent = (1e6) * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec);

    printf("%d ms\n",time_spent);    
    return 0;
}
