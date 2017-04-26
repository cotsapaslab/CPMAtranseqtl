
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <time.h>

#define MIN(a,b) ((a)<(b)?(a):(b))

char *gettime() {
  time_t now;
  time(&now);
  return ctime(&now);
}

int main(int argc, char *argv[])
{
  int ifd;
  FILE *ifp, *ofp;
  long rows, cols;
  long window=atol(argv[3]);
  long startcol=atoi(argv[4]);

  float *map, *buf, *tmpptr;
  int c, r, cc, did;
  float val;
  float *ptrs[window];
  long i, thiswin, doing, dxid, expected, got;

  ifp = fopen(argv[1], "rb");
  fseek(ifp, -2*sizeof(int), SEEK_END);
  fread(&rows, sizeof(int), 1, ifp);
  fread(&cols, sizeof(int), 1, ifp);

  got=ftell(ifp);
  expected=(rows*cols*sizeof(float)+2*sizeof(int));
  printf("rows %d cols %d\n", rows, cols);
  if (got != expected) {
    printf("Size issue: expected %d got %d", expected, got);
    exit(-1);
  }

  fclose(ifp);
  if ((access(argv[2], W_OK)) != -1)
    ofp=fopen(argv[2], "r+b");
  else
    ofp=fopen(argv[2], "wb");

  fseek(ofp, startcol*rows*sizeof(float), SEEK_SET);

  ifd=open(argv[1], O_RDONLY);
  
  buf = calloc(rows*window, sizeof(float));
  if (!buf) {
    printf("calloc failed\n");
    exit(-1);
  }

  map=mmap(0, rows*cols*4, PROT_READ, MAP_PRIVATE, ifd, 0);
  if (map==MAP_FAILED) {
    printf("Mmap failed\n");
    exit(-1);
  }

  for (c=startcol; c<cols; c+=window) {
    printf("doing col %d %s\n", c, gettime());
    // initialize ptrs into buffer
    for (i=0, tmpptr=buf; i<window; ++i, tmpptr+=rows) 
      ptrs[i]=tmpptr;

    thiswin=MIN(window, cols-c);
    for (r=0; r<rows; ++r) {
      for (cc=0; cc<thiswin; ++cc) 
	*ptrs[cc]++ = map[r*cols+c+cc];
    }

    doing=sizeof(float)*rows*thiswin;
    did=fwrite(buf, doing, 1, ofp);
    if (did!=1) {
      printf("fwrite failed loser %d\n", did);
      exit(-1);
    }
  }
  // intentionally reversed, because in the new file, rows and cols are reversed.
  fwrite(&cols, sizeof(int), 1, ofp);
  fwrite(&rows, sizeof(int), 1, ofp);
  exit(0);
}
