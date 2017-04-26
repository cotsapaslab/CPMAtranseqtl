
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <time.h>

int main(int argc, char *argv[])
{
  int ifd;
  FILE *ifp;
  int rows, cols;

  long at=atol(argv[2]);

  ifp = fopen(argv[1], "r+b");
  fseek(ifp, at, SEEK_SET);

  fread(&rows, sizeof(int), 1, ifp);
  fread(&cols, sizeof(int), 1, ifp);

  printf("at %ld rows %d cols %d", at, rows, cols);
  exit(0);
}
