#include "getSequenceCorrelations.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hashBases.hpp"

char chrNames[23][3];
float valueTable[23][4096],freqValueTable[3][4096], corrs[23][23];
int totals[4096][2],vtotals[23][3];
char baseTable[4096][2][7], backgroundBaseTable[1024][2][6];

baseHasher hasher;

double correl(float* x, float* y, int count)
{
	double xbar, ybar, bxy, byx, sigmaxy = 0, sigmax2 = 0, sigmay2 = 0, sigmax = 0, sigmay = 0, r;
	int i;
	for (i = 0; i < count; ++i)
	{
		sigmax += x[i];
		sigmay += y[i];
		sigmaxy += x[i] * y[i];
	}
	xbar = sigmax / count;
	ybar = sigmay / count;
	byx = bxy = 0;
	for (i = 0; i < count; ++i)
	{
		bxy += (x[i] - xbar) * (y[i] - ybar);
		sigmax2 += x[i] * x[i];
		sigmay2 += y[i] * y[i];
	}
	byx = bxy;
	bxy /= sigmax2 - sigmax * sigmax / count;
	byx /= sigmay2 - sigmay * sigmay / count;
	r = sqrt(bxy * byx) * (bxy < 0 ? -1 : 1);
	return r;
}

int getCorrs(int count,int nCols=23) {
	int cx, cy;
	for (cx = 0; cx < nCols; ++cx) {
		corrs[cx][cx] = 1.0;
		for (cy = cx + 1; cy < nCols; ++cy) {
			corrs[cx][cy] = corrs[cy][cx] = correl(valueTable[cx], valueTable[cy], count);
		}
	}
	return 1;
}

int printCorrs(FILE* fo,int nCols=23)
{
	fprintf(fo, nCols==23?"CHR":"Column");
	int cx, cy;
	for (cx = 0; cx < nCols; ++cx)
		fprintf(fo, "\t%s",chrNames[cx]);
	fprintf(fo, "\n");
	for (cy = 0; cy < nCols; ++cy) {
		if (nCols==23)
			fprintf(fo, "%s", chrNames[cy]);
		else
			fprintf(fo, "%d", cy);
		for (cx = 0; cx < nCols; ++cx)
			fprintf(fo, "\t%.3f", corrs[cx][cy]);
		fprintf(fo, "\n");
	}
	fprintf(fo, "\n");
	return 1;
}

int main(int argc, char* argv[])
{
	int i, c, f, l, b;
	unsigned long sum; // signed long is not big enough!!
	float c1, c2,v1,v2;
	FILE* fi, * fo, * ft;
	char* chr, line[1000], fn[100], seq[7], compSeq[7],countFileRoot[100],
		backgroundCountsFileRoot[100], compPairCountsFileRoot[100],
		compPairBackgroundCountsFileRoot[100], compPairCorrectedCountsFileRoot[100],startingBase;
	unsigned int hash;
	for (c = 0; c < 22; ++c)
		sprintf(chrNames[c], "%d", c + 1);
	strcpy(chrNames[22], "X");
	sprintf(fn, "sequenceCorrelations.txt");
	fo = fopen(fn, "w");
	for (i = 0; i < 4096; ++i)
		totals[i][0] = totals[i][1] = 0;
	for (c = 0; c < 23; ++c) {
		chr = chrNames[c];
		sprintf(backgroundCountsFileRoot, "backgroundCounts.%s", chr);
		sprintf(fn, "%s.txt", backgroundCountsFileRoot);
		fi = fopen(fn, "r");
		i = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %f", baseTable[i][0], &valueTable[c][i]);
			if (c!=18)
				totals[i][0] += valueTable[c][i];
			totals[i][1] += valueTable[c][i];
			++i;
		}
		fclose(fi);
	}
	ft = fopen("backgroundCounts.not19.txt", "w");
	for (l = 0; l < i; ++l)
		fprintf(ft, "%s\t%d\n", baseTable[l][0], totals[l][0]);
	fclose(ft);
	ft = fopen("backgroundCounts.total.txt", "w");
	sum = 0;
	for (l = 0; l < i; ++l) {
		fprintf(ft, "%s\t%d\n", baseTable[l][0], totals[l][1]);
		sum += totals[l][1];
	}
	fclose(ft);
	fprintf(fo, "Total number of background sequences = %lu\n\n", sum);
	getCorrs(i);
	fprintf(fo, "Correlations for background counts\n\n");
	printCorrs(fo);
	for (f = 0; f < 3; ++f) {
		for (i = 0; i < 4096; ++i)
			totals[i][0]= 0;
		for (c = 0; c < 23; ++c) {
			chr = chrNames[c];
			sprintf(countFileRoot, "counts.%s.%d", chr, f);
			sprintf(fn, "%s.txt", countFileRoot);
			fi = fopen(fn, "r");
			i = 0;
			while (fgets(line, 999, fi)) {
				sscanf(line, "%s %f", baseTable[i][0], &valueTable[c][i]);
				freqValueTable[f][i] += valueTable[c][i];
				vtotals[c][f] += valueTable[c][i];
				totals[i][0]+= valueTable[c][i];
				++i;
			}
			fclose(fi);
			sprintf(fn, "counts.total.%d.txt", f);
			ft = fopen(fn, "w");
			sum = 0;
			for (l = 0; l < i; ++l) {
				fprintf(ft, "%s\t%d\n", baseTable[l][0], totals[l][0]);
				sum += totals[l][0];
			}
			fclose(ft);
		}
		fprintf(fo, "Total number of %d variants = %lu\n\n", f, sum);
		getCorrs(i);
		fprintf(fo, "Correlations for variant counts %d\n\n", f);
		printCorrs(fo);
	}
	ft = fopen("variantTotalsByFrequency.txt", "w");
	for (c = 0; c < 23; ++c) {
		fprintf(ft, "%s", chrNames[c]);
		for (f = 0; f < 3; ++f)
			fprintf(ft, "\t%d", vtotals[c][f]);
		fprintf(ft, "\n");
	}
	fclose(ft);
	for (l = 0; l < i; ++l)
		for (f = 0; f < 3; ++f)
			valueTable[f][l] = freqValueTable[f][l];
	getCorrs(i, 3);
	fprintf(fo, "Correlations between variant frequency categories\n\n");
	printCorrs(fo, 3);

	for (f = 0; f < 3; ++f) {
		for (i = 0; i < 4096; ++i)
			totals[i][0] = 0;
		for (c = 0; c < 23; ++c) {
			chr = chrNames[c];
			sprintf(countFileRoot, "correctedCounts.%s.%d", chr, f);
			sprintf(fn, "%s.txt", countFileRoot);
			fi = fopen(fn, "r");
			i = 0;
			while (fgets(line, 999, fi)) {
				sscanf(line, "%s %f", baseTable[i][0], &valueTable[c][i]);
				totals[i][0] += valueTable[c][i];
				++i;
			}
			fclose(fi);
		}
		getCorrs(i);
		fprintf(fo, "Correlations for corrected variant counts %d\n\n", f);
		printCorrs(fo);
		sprintf(fn, "correctedCounts.total.%d.txt",f);
		ft = fopen(fn, "w");
		for (l = 0; l < i; ++l)
			fprintf(ft, "%s\t%d\n", baseTable[l][0], totals[l][0]);
		fclose(ft);
	}
	
	ft = fopen("correctedCounts.total.byFrequency.txt", "w");
	for (l = 0; l < i; ++l)
		fprintf(ft, "%s\t%d\t%d\t%d\n", baseTable[l][0], vtotals[l][0], vtotals[l][2], vtotals[l][2]);
	fclose(ft);
	for (i = 0; i < 4096; ++i)
		totals[i][0] = totals[i][1] = 0;
	for (c = 0; c < 23; ++c) {
		chr = chrNames[c];
		sprintf(compPairBackgroundCountsFileRoot, "compPairBackgroundCounts.%s", chr);
		sprintf(fn, "%s.txt", compPairBackgroundCountsFileRoot);
		fi = fopen(fn, "r");
		i = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %s %f %f", backgroundBaseTable[i][0], backgroundBaseTable[i][1], &c1, &c2);
			valueTable[c][i] = log(c1 / c2);
			totals[i][0] += c1;
			totals[i][1] += c2;
			++i;
		}
		fclose(fi);
	}
	getCorrs(i);
	fprintf(fo, "Correlations for background count strand asymmetries\n\n");
	printCorrs(fo);
	ft = fopen("compPairBackgroundCounts.total.txt", "w");
	for (l = 0; l < i; ++l)
		fprintf(ft, "%s\t%s\t%d\t%d\n", backgroundBaseTable[l][0], backgroundBaseTable[l][1], totals[l][0], totals[l][1]);
	fclose(ft);
	for (f = 0; f < 3; ++f) {
		for (i = 0; i < 4096; ++i)
			totals[i][0] = totals[i][1] = 0;
		for (c = 0; c < 23; ++c) {
			chr = chrNames[c];
			sprintf(compPairCorrectedCountsFileRoot, "compPairCorrectedCounts.%s.%d", chr,f);
			sprintf(fn, "%s.txt", compPairCorrectedCountsFileRoot);
			fi = fopen(fn, "r");
			i = 0;
			while (fgets(line, 999, fi)) {
				sscanf(line, "%s %s %f %f", baseTable[i][0], baseTable[i][1], &c1, &c2);
				valueTable[c][i] = log(c1 / c2);
				totals[i][0] += c1;
				totals[i][1] += c2;
				++i;
			}
			fclose(fi);
		}
		getCorrs(i);
		fprintf(fo, "Correlations for corrected variant count strand asymmetries %d\n\n",f);
		printCorrs(fo);
		sprintf(fn, "compPairCorrectedCounts.total.%d", f);
		ft = fopen(fn, "w");
		for (l = 0; l < i; ++l)
			fprintf(ft,"%s\t%s\t%d\t%d\n", baseTable[l][0], baseTable[l][1], totals[l][0], totals[l][1]);
		fclose(ft);
	}
	for (f = 0; f < 3; ++f) {
		for (i = 0; i < 4096; ++i)
			totals[i][0] = totals[i][1] = 0;
		sprintf(fn, "compPairCorrectedCounts.total.%d", f);
		fi = fopen(fn, "r");
		i = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %s %f %f", baseTable[i][0], baseTable[i][1], &c1, &c2);
			hash = hasher.hashBases(baseTable[i][0], 6);
			freqValueTable[0][hash] = c1;
			freqValueTable[1][hash] = c2;
			++i;
		}
		fclose(fi);
		fi = fopen("compPairBackgroundCounts.total.txt", "r");
		i = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %s %f %f", backgroundBaseTable[i][0], backgroundBaseTable[i][1], &c1, &c2);
			v1 = v2 = 0;
			sprintf(seq, "%s%c", backgroundBaseTable[i][0], backgroundBaseTable[i][0][2]);
			// for this bit you will need to see whether to take the complement of the variant
			for (b = 0; b < 4; ++b) {
				seq[2] = baseNames[b];
				if (seq[2] == seq[5])
					continue;
				if (seq[2] == 'C' || seq[2] == 'T') {
					hash = hasher.hashBases(seq, 6);
					v1 += freqValueTable[0][hash];
					v2 += freqValueTable[1][hash];
				}
				else {
					hasher.getComplement(compSeq, seq, 5);
					hasher.getComplement(compSeq+5, seq+5, 1);
					hash = hasher.hashBases(compSeq,6);
					v2 += freqValueTable[0][hash];
					v1 += freqValueTable[1][hash];
				}
			}
			valueTable[0][i] = log(c1 / c2);
			valueTable[1][i] = log(v1 / v2);
			++i;
		}
		fclose(fi);
		fprintf(fo, "Correlation between strand mismatch for background sequences and for their three generating %d variants: R = %.4f\n", f,correl(valueTable[0], valueTable[1], i));
	}
	fclose(fo);

	return 0;
}