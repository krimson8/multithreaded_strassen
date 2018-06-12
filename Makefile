all: strassen.c
	gcc -Wall -o main strassen.c -lpthread

clean:
	rm main
