cosamp.exe: cosamp.c
	gcc -O3 -o cosamp cosamp.c f2c.h clapack.h