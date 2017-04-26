
all: transpose2 fixCounts getCounts

transpose2: transpose2.c
	$(CC) -o transpose2 transpose2.c

getCounts: getCounts.c
	$(CC) -o getCounts getCounts.c

fixCounts: fixCounts.c
	$(CC) -o fixCounts fixCounts.c


