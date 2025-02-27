all: mat3 brin

mat3: mat3.c
	gcc -g -o  mat3 mat3.c

brin: brin.c
	gcc -g -o  brin brin.c
	
clean:
	rm -f *.o *~ a.out mat3 brin
