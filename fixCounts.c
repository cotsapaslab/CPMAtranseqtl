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

  long at=atol(argv[2]);
  long rows=atoi(argv[3]);
  long cols=atoi(argv[4]);

  long expected, got;

  ifp = fopen(argv[1], "r+b");
  fseek(ifp, at, SEEK_SET);
  got=ftell(ifp);

  expected=(rows*cols*sizeof(float));
  if (got != expected) {
    printf("Size issue: expected %d got %d", expected, got);
    exit(-1);
  }

  fwrite(&rows, sizeof(int), 1, ifp);
  fwrite(&cols, sizeof(int), 1, ifp);

  exit(0);
}
