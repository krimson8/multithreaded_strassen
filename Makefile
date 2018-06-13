all: strassen.c
	gcc -Wall -o main strassen.c -lpthread

opt: strassen.c
	gcc -Wall -O3 -o main strassen.c -lpthread

clean:
	rm main
