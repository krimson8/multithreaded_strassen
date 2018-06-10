#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <iostream>

using namespace std;

int n;
int ma[4096][4096];
int mb[4096][4096];
long int ans[4096][4096];

int thread_count = 4;
pthread_mutex_t mutexA = PTHREAD_MUTEX_INITIALIZER;
int counter = 0;

void* child(void *pos){
    cout << "Thread number: " << pthread_self () << "\n";
    int m = n / 2;

    int p = *(int*)pos;
    int px = p / 2;
    int py = p % 2;
    for(int i = px * m; i < m + px * m; i++) {
        for(int j = py * m; j < m + py * m; j++) {
            ans[i][j] = 0;
            for(int k = 0; k < n; k++) {
                ans[i][j] += ma[i][k] * mb[k][j];
            }
        }
    }

    pthread_mutex_lock (&mutexA);
    counter++;
    pthread_mutex_unlock (&mutexA);

    pthread_exit(NULL);
}

int main() {
    FILE *fp = fopen("test1", "r");

    // read matrix A
    fscanf(fp, "%d %d", &n, &n);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            fscanf(fp, "%d", &ma[i][j]);
        }
    }

    // read matrix B
    fscanf(fp, "%d %d", &n, &n);
    for(int i=0 ; i < n; i++) {
        for(int j = 0; j < n; j++) {
            fscanf(fp, "%d", &mb[i][j]);
        }
    }
    fclose(fp);

    clock_t begin = clock();

    pthread_t t[thread_count];
    for(int i = 0; i < thread_count; i++) {
        pthread_create(&t[i], NULL, child, &i);
    }

    for(int i = 0; i < thread_count; i++) {
        pthread_join(t[i], NULL);
    }
    clock_t end = clock();

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%ld ", ans[i][j]);
        }
        printf("\n");
    }

    double time_spent = (double)(end - begin) / (CLOCKS_PER_SEC*4);
    printf("%lf\n",time_spent);    
}