OBJS = GetBaseCounts.o VariantFile.o VariantEntry.o Global.o
CC = g++
CFLAGS = -Wall -c
LFLAGS = -Wall -o3
GZFLAGS = -I./gzstream/include/ -L./gzstream/lib/ -lgzstream -lz
BAMFLAGS = -I./bamtools-master/include/ -L./bamtools-master/lib/ -lbamtools
OMPFLAGS = -fopenmp

GetBaseCounts: $(OBJS)
	$(CC) $(OBJS) $(LFLAGS) -o GetBaseCounts $(GZFLAGS) $(BAMFLAGS) $(OMPFLAGS)

GetBaseCounts.o: GetBaseCounts.cpp VariantFile.h VariantEntry.h Global.h
	$(CC) $(CFLAGS) $(GZFLAGS) $(BAMFLAGS) $(OMPFLAGS) GetBaseCounts.cpp

VariantFile.o: VariantFile.cpp VariantFile.h
	$(CC) $(CFLAGS) $(GZFLAGS) VariantFile.cpp

VariantEntry.o: VariantEntry.cpp VariantEntry.h
	$(CC) $(CFLAGS) VariantEntry.cpp

Global.o: Global.cpp Global.h
	$(CC) $(CFLAGS) $(OMPFLAGS) Global.cpp

clean:
	rm -f *.o GetBaseCounts

