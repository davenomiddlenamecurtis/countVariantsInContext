#include "countVariantsInContext.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

baseHasher hasher;

const char* chr = "22";
int countTable[3][4096];

char outputSequence[7];

int outputCounts(FILE* fo, int level,int *counts) {
	int i;
	if (level == 0) {
		for (i = 0; i < 4; ++i) {
			outputSequence[5 - level] = baseNames[i];
			fprintf(fo, "%s\t%d\n", outputSequence, counts[hasher.hashBases(outputSequence, 6)]);
		}
		return 1;
	}
	else {
		for (i = 0; i < 4; ++i) {
			outputSequence[5 - level] = baseNames[i];
			outputCounts(fo, level - 1,counts);
		}
		return 1;
	}
	return 0;
}


int main(int argc, char* argv[])
{
	int i,start,AC;
	unsigned int hash;
	FILE* fi,*fo;
	faSequenceFile fa;
	char inFile[100],outFileRoot[100],refFile[100];
	if (argc > 1) {
		chr = argv[1];
	}
	sprintf(inFile, "ukb24308_c%s_b0_v1.head.pvar", chr);
	sprintf(outFileRoot, "counts.%s", chr);
	sprintf(refFile, "CHR%s.FA", chr);
	fa.init(refFile);
	fi = fopen(inFile, "r");
	char* posPtr, * ptr, line[1000],sequence[6],fn[100];
	fgets(line, 999, fi);
	while (fgets(line, 999, fi)) {
		ptr = line;
		while (!isspace(*ptr++));
		while (isspace(*ptr++));
		posPtr = ptr-1;
		for (i = 0; i < 3; ++i)
			while (*ptr++ != ':');
		if (ptr[1]!=':' || !isspace(ptr[3]))
			continue; // only use SNVs
		start = atoi(posPtr) - 2;
		fa.getSequence(sequence, start, 5);
		if (sequence[0] == 'N')
			continue;
		sequence[5] = ptr[2];
		hash=hasher.hashBases(sequence,6);
		ptr += 7;
		while (*ptr++ != 'A' || *ptr++ != 'C') // look for AC= string in line
			;
		AC = atoi(ptr + 1);
		if (AC>3)
			++countTable[0][hash];
		else if (AC >1)
			++countTable[1][hash];
		else if (AC ==1)
				++countTable[2][hash];
// 22 10510001 DRAGEN:chr22:10510001:G:T G T 0 LowGTR AC=2;AN=17868;AF=0.000111932;NS=490541;NS_GT=8934;NS_NOGT=7484;NS_NODATA=474123
	}
	fclose(fi);
	for (i = 0; i < 3; ++i)
	{
		sprintf(fn, "%s.%d.txt", outFileRoot, i);
		fo = fopen(fn, "w");
		outputCounts(fo, 5,countTable[i]);
		fclose(fo);
	}

	printf("%d\n", hasher.hashBases("TGCA", 4));
	printf("%d\n", hasher.hashBases("AAAA", 4));
	printf("%d\n", hasher.hashBases("AAAT", 4));
	printf("%d\n", hasher.hashBases("TAAA", 4));
	printf("%d\n", hasher.hashBases("TTTT", 4));
	printf("%d\n", hasher.hashBases("CCCC", 4));
	printf("%d\n", hasher.hashBases("GGGG", 4));
	return 0;
}