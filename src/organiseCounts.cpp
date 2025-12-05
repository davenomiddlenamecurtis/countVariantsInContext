#include "organiseCounts.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

baseHasher hasher;

const char* chr = "22";
int countTable[3][4096], backgroundCountTable[1024],correctedCountTable[3][4096];
char baseTable[4096][7],backgroundBaseTable[1024][6];
float correction[1024];

int main(int argc, char* argv[])
{
	int i,c,count,diffHash[2],tot,correctedDiffHash[2];
	float maxDiff, diff, maxDiffCorrected;
	unsigned int hash,background;
	FILE* fi,*fo,*fc;
	char line[1000],fn[100],seq[7], countFileRoot[100], 
		backgroundCountsFileRoot[100], compPairCountsFileRoot[100],correctedCountsFileRoot[100],
		compPairBackgroundCountsFileRoot[100], compPairCorrectedCountsFileRoot[100];
	if (argc > 1) {
		chr = argv[1];
	}
	sprintf(countFileRoot, "counts.%s", chr);
	sprintf(backgroundCountsFileRoot, "backgroundCounts.%s", chr);
	sprintf(compPairBackgroundCountsFileRoot, "compPairBackgroundCounts.%s", chr);
	sprintf(compPairCountsFileRoot, "compPairCounts.%s", chr);
	sprintf(correctedCountsFileRoot, "correctedCounts.%s", chr);
	sprintf(compPairCorrectedCountsFileRoot, "compPairCorrectedCounts.%s", chr);
	sprintf(fn, "%s.txt", backgroundCountsFileRoot);
	fi = fopen(fn, "r");
	seq[5] = seq[6] = '\0';
	while (fgets(line, 999, fi)) {
		sscanf(line, "%s %d", seq, &count);
		hash = hasher.hashBases(seq, 5);
		backgroundCountTable[hash] = count;
		strcpy(backgroundBaseTable[hash], seq);
	}
	fclose(fi);
	tot = 0;
	for (i = 0; i < 1024; ++i)
		tot += backgroundCountTable[i];
	for (i = 0; i < 1024; ++i)
		correction[i]=(tot/1024.0)/backgroundCountTable[i];
	sprintf(fn, "%s.txt", compPairBackgroundCountsFileRoot);
	fo = fopen(fn, "w");
	maxDiff = 0;
	seq[5] = seq[6] = '\0';
	for (i = 0; i < 1024; ++i) {
		if (backgroundBaseTable[i][2] == 'C' || backgroundBaseTable[i][2] == 'T') {
			hasher.getComplement(seq, backgroundBaseTable[i], 5);
			hash = hasher.hashBases(seq, 5);
			diff = fabs(log(backgroundCountTable[i]/(float)backgroundCountTable[hash]));
			if (diff > maxDiff) {
				maxDiff = diff;
				diffHash[0] = i;
				diffHash[1] = hash;
			}
			fprintf(fo, "%s\t%s\t%d\t%d\n", backgroundBaseTable[i], seq, backgroundCountTable[i], backgroundCountTable[hash]);
		}
	}
	fclose(fo);
	printf("Maximum ratio for counts between complementary pairs for background sequences is between %d and %d for %s and %s\n",
		backgroundCountTable[diffHash[0]], backgroundCountTable[diffHash[1]],
		backgroundBaseTable[diffHash[0]],backgroundBaseTable[diffHash[1]]);
	// This paper seems to say that strand asymmetries can occur at the chromosomal level
	// https://pmc.ncbi.nlm.nih.gov/articles/PMC10030826/

	// means best to keep sequences separate till after correction for background frequency
	// so one can see strand-specific effects on mutations

	for (c = 0; c < 3; ++c) {
		sprintf(fn, "%s.%d.txt", countFileRoot, c);
		fi = fopen(fn, "r");
		sprintf(fn, "%s.%d.txt", correctedCountsFileRoot, c);
		fo = fopen(fn, "w");
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %d", seq, &count);
			hash = hasher.hashBases(seq, 6);
			countTable[c][hash] = count;
			if (c == 0) {
				strcpy(baseTable[hash], seq);
			}
			background = hash >> 2;
			correctedCountTable[c][hash] = countTable[c][hash] * correction[background];
			if (seq[2]!=seq[5])
				fprintf(fo, "%s\t%d\n", seq, correctedCountTable[c][hash]);
		}
		fclose(fi);
		fclose(fo);
	}
	for (c = 0; c < 3; ++c) {
		sprintf(fn, "%s.%d.txt", compPairCountsFileRoot, c);
		fo = fopen(fn, "w");
		sprintf(fn, "%s.%d.txt", compPairCorrectedCountsFileRoot, c);
		fc = fopen(fn, "w");
		maxDiff = maxDiffCorrected= 0;
		seq[6] = '\0';
		for (i = 0; i < 4096; ++i) {
			if ((baseTable[i][2] == 'C' || baseTable[i][2] == 'T') && baseTable[i][2]!= baseTable[i][5]) {
				hasher.getComplement(seq, baseTable[i], 5);
				hasher.getComplement(seq+5, baseTable[i]+5, 1);
				hash = hasher.hashBases(seq, 6);
				diff = fabs(log(countTable[c][i]/(float)countTable[c][hash]));
				if (diff > maxDiff) {
					maxDiff = diff;
					diffHash[0] = i;
					diffHash[1] = hash;
				}
				diff = fabs(log(correctedCountTable[c][i]/(float)correctedCountTable[c][hash]));
				if (diff > maxDiffCorrected) {
					maxDiffCorrected = diff;
					correctedDiffHash[0] = i;
					correctedDiffHash[1] = hash;
				}
				fprintf(fo, "%s\t%s\t%d\t%d\n", baseTable[i], seq, countTable[c][i], countTable[c][hash]);
				fprintf(fc, "%s\t%s\t%d\t%d\n", baseTable[i], seq, correctedCountTable[c][i], correctedCountTable[c][hash]);
			}
		}
			fclose(fo);
			fclose(fc);
			printf("Maximum ratio for %d counts between complementary pairs for sequences is between %d and %d for %s and %s\n",
				c,
				countTable[c][diffHash[0]], countTable[c][diffHash[1]],
				baseTable[diffHash[0]], baseTable[diffHash[1]]);
			printf("Maximum ratio for %d corrected counts between complementary pairs for sequences is between %d and %d for %s and %s\n",
				c,
				correctedCountTable[c][correctedDiffHash[0]], correctedCountTable[c][correctedDiffHash[1]],
				baseTable[correctedDiffHash[0]], baseTable[correctedDiffHash[1]]);
	}
	return 0;
}