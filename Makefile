CC      = clang
CFLAGS  = -fopenmp        # compiler flags (for compilation)
LDFLAGS = -fopenmp        # linker flags (for linking)
OBJ     = main.o guassonFilter.o blur.o
EXEC    = main

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(LDFLAGS) -o $(EXEC) $(OBJ)

main.o: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o

guassonFilter.o: guassonFilter.c
	$(CC) $(CFLAGS) -c guassonFilter.c -o guassonFilter.o

blur.o: blur.c
	$(CC) $(CFLAGS) -c blur.c -o blur.o

clean:
	rm -f $(EXEC) $(OBJ)
