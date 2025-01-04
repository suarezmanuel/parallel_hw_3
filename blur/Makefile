CC      = gcc
CFLAGS  = -g -fopenmp -msse -msse2 -march=native -O2        # compiler flags (for compilation)
LDFLAGS = -fopenmp -lm        # linker flags (for linking)
OBJ     = main.o guassonFilter.o blur.o
EXEC    = main

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(LDFLAGS) -o $(EXEC) $(OBJ)
	rm -f blur.o guassonFilter.o main.o

main.o: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o

guassonFilter.o: guassonFilter.c
	$(CC) $(CFLAGS) -c guassonFilter.c -o guassonFilter.o

blur.o: blur.c
	$(CC) $(CFLAGS) -c blur.c -o blur.o

run:
	./main

clean:
	rm -f $(EXEC) $(OBJ)
