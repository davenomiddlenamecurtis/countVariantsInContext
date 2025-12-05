#include "getBackgroundCounts.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

baseHasher hasher;

int backgroundCountTable[1024];
const char* chr = "22";


// background sequence consists of 5 bases

char outputSequence[6];

int outputBackgroundCounts(FILE* fo, int level,int *counts) {
	int i;
	if (level == 0) {
		for (i = 0; i < 4; ++i) {
			outputSequence[4 - level] = baseNames[i];
			fprintf(fo, "%s\t%d\n", outputSequence, counts[hasher.hashBases(outputSequence, 5)]);
		}
		return 1;
	}
	else {
		for (i = 0; i < 4; ++i) {
			outputSequence[4 - level] = baseNames[i];
			outputBackgroundCounts(fo, level - 1,counts);
		}
		return 1;
	}
	return 0;
}

int main(int argc, char* argv[])
{
	int i,c,tot;
	unsigned int hash;
	char fn[100];
	const char *ptr;
	FILE* fi,*fo;
	if (argc > 1) {
		chr = argv[1];
	}
	sprintf(fn, "CHR%s.FA", chr);
	fi = fopen(fn, "r");
	i = 0;
	while ((c = fgetc(fi))!=EOF) {
		if (c == 'N' || isspace(c) || !strchr(baseNames, c))
			continue;
		outputSequence[i++] = c;
		if (i > 4)
			break;
	}
	// now have first lot of five letters
	hash = hasher.hashBases(outputSequence, 5);
	++backgroundCountTable[hash];
	tot = 1;
	while ((c = fgetc(fi))!=EOF) {
		if (c == 'N' || isspace(c) || (ptr=strchr(baseNames, c))==0)
			continue;
		hash &= 255;
		hash <<= 2;
		hash |= (ptr - baseNames);
		++backgroundCountTable[hash];
		++tot;
	}

	fclose(fi);
	sprintf(fn, "backgroundCounts.%s.txt", chr);
	fo = fopen(fn, "w");
	outputBackgroundCounts(fo, 4,backgroundCountTable);
	fclose(fo);
	printf("Dealt with %d bases\n", tot);
	return 0;
}