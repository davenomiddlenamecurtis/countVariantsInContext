#include "modelStrandBias.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "hashBases.hpp"
#include "runModels.hpp"

#if 0
in general, each row will have a sequence, a count and a probability
when I analyse corrected counts, the total count will be the number of times background occurred
(same for every sequence)
the names of component sequences will be generated, then betas will be 0 or 1 depending on whether they match full sequence
each sequence represents a complementary pair
#endif

#define MAXLEVELS 10 // for pseudo-recursion
// level will mean the opposite of what it did in outputCounts() - it means closer to the end
int generateSequence(char* s, int len,int start)
{
	static int level, i[MAXLEVELS],j;
	static char bases[MAXLEVELS + 1];
	if (start) {
		for (j = 0; j < len; ++j) {
			bases[j] = baseNames[0];
			i[j] = 0;
		}
		i[j-1] = -1;
		return 1;
	}
	else {
		if (i[len-1] < 3) {
			++i[len - 1];
			bases[len - 1] = baseNames[i[len - 1]];
			strcpy(s, bases);
			return 1;
		}
		else {
			for (j = len - 1; j > 0; --j)
				if (i[j] < 3)
					break;
			if (j == 0 && i[j] == 3)
				return 0;
			else {
				++i[j];
				bases[j] = baseNames[i[j]];
				for (++j; j < len; ++j) {
					i[j] = 0;
					bases[j] = baseNames[i[j]];
				}
				strcpy(s, bases);
				return 1;
			}
		}
	}
}

// just add up counts to get averages for each beta
void setStartingBetasFromCounts(glfModel* m)
{
	float N=0, *sigmaX,p;
	int b, r;
	assert((sigmaX = (float*)calloc(m->nCol, sizeof(float))) != 0);
	for (r = 0; r < m->nRow; ++r) {
		N += m->F[r];
		for (b = 0; b < m->nCol; ++b) 
			if (m->toFit[b]) {
				sigmaX[b] += m->X[r][b] * m->F[r];
			}
	}
	for (b = 0; b < m->nCol; ++b)
		if (m->toFit[b]) {
			p=sigmaX[b]/N;
			m->beta[b] = log(p / (1 - p));
			m->SE[b] = log(sqrt(N * p * (1 - p))); 
			//actually log of SE of p, not odds, but I do not care
		}
	free(sigmaX);
}

// produce variants in same order as is used by SigProfilerAssignment for input matrix (FFS)
int generateSBSVariantForInputMatrix(char* s, int start)
{
	static int i, j, k, l;
	static const char* CT = "CT", * target[2] = { "AGT","ACG" }, * ACGT = "ACGT";
	if (start == 1) {
		i = j = k = 0;
		l = -1;
		return 1;
	}

	++l;
	if (l > 3) {
		l = 0;
		++k;
		if (k <= 3 && ACGT[k] == CT[j])
			++k;
		if (k > 3) {
			k = 0;
			++j;
			if (j > 1) { // C or T
				j = 0;
				++i;
				if (i > 3)
					return 0;
			}
		}
	}
	sprintf(s, "_%c%c%c_%c", ACGT[i], CT[j], ACGT[l], ACGT[k]);
	return 1;
}

// produce variants in same order as is output by SigProfilerAssignment (FFS)
int generateSBSVariant(char* s, int start)
{
	static int i, j, k, l;
	static const char* CT = "CT", * target[2] = { "AGT","ACG" }, * ACGT = "ACGT";
	if (start == 1) {
		i = j = k = 0;
		l = -1;
		return 1;
	}

	++l;
	if (l > 3) {
		l = 0;
		++k;
		if (k > 3) {
			k = 0;
			++j;
			if (j > 2) {
				j = 0;
				++i;
				if (i > 1)
					return 0;
			}
		}
	}
	sprintf(s, "_%c%c%c_%c", ACGT[k], CT[i], ACGT[l], target[i][j]);
	return 1;
}

int matches(char* seq, char* comp, int len)
{
	int i;
	for (i = 0; i < len; ++i)
		if (comp[i] != '_' && comp[i] != seq[i])
			return 0;
	return 1;
}

int getBetas(double * beta, glfModel* m)
{
	double incBeta0=0;
	int b;
	if (m->isNormalised)
	{
		incBeta0 = 0;
		for (b = 0; b < m->nCol; ++b)
		{
			if (m->SD[b] && m->toUse[b])
			{
				incBeta0 -= m->beta[b] * m->mean[b] / m->SD[b];
				beta[b] = m->beta[b] / m->SD[b];
			}
			else
				beta[b] = 0;
		}
		beta[m->nCol] = m->beta[m->nCol] + incBeta0;
	}
	else
		for (b = 0; b < m->nCol + 1; ++b) {
			if (m->toUse[b])
				beta[b] = m->beta[b];
			else
				beta[b] = 0;
		}
	return 1;
}

void printModel(FILE* fo, const char* LLstr, double LL, glfModel* m)
{
	// change this to be row-wise and use names
	int b, bb, c;
	double* realBeta, * realSE;
	double incBeta0, incSEBeta0;
	realBeta = (double*)calloc(m->nCol + 1, sizeof(double));
	realSE = (double*)calloc(m->nCol + 1, sizeof(double));
	fprintf(fo, "%s = %.2f\n", LLstr, LL);
	fprintf(fo, "beta\tvalue\tSE\tz\n");
	if (m->isNormalised)
	{
		incSEBeta0 = incBeta0 = 0;
		for (c = 0; c < m->nCol; ++c)
		{
			if (m->SD[c] && m->toUse[c])
			{
				incBeta0 -= m->beta[c] * m->mean[c] / m->SD[c];
				incSEBeta0 -= m->SE[c] * m->mean[c] / m->SD[c];
				realBeta[c] = m->beta[c] / m->SD[c];
				realSE[c] = m->SE[c] / m->SD[c];
			}
		}
		realBeta[m->nCol] = m->beta[m->nCol] + incBeta0;
		realSE[m->nCol] = m->SE[m->nCol] + incSEBeta0;
	}
	else
		for (b = 0; b < m->nCol + 1; ++b)
		{
			realBeta[b] = m->beta[b];
			realSE[b] = m->SE[b];
		}
	for (b = 0; b < m->nCol + 1; ++b)
	{
		bb = (b + m->nCol) % (m->nCol + 1); // print last first
		if (m->toUse[bb])
			fprintf(fo, "%s\t%.5f\t%.5f\t%.5f\n", m->name[bb], realBeta[bb], realSE[bb], m->toFit[bb] ? realBeta[bb] / realSE[bb] : 0.0);
	}
	fprintf(fo, "\n");
	free(realBeta);
	free(realSE);
}

float evaluateModel(FILE* fo, glfModel* m, const char* name,int *useThese,double *startingBetas,int useLinearRegression)
{
	float lnL;
	int b;
	for (b = 0; b < m->nCol+1; ++b)
	{
		m->toUse[b] = useThese?useThese[b]:1;
		m->beta[b] = startingBetas? startingBetas[b]:0;
		m->toFit[b] = useThese ? useThese[b] : 1;
	}
	m->normalise();
	if (useLinearRegression)
		m->useLinearRegression(1);
	lnL = m->maximiseLnL();
	m->getSEs();
	//	m->deNormalise(); // avoid repeated normalising and deNormalising
	if (fo)
		printModel(fo, name, lnL, m);
	return lnL;
}


