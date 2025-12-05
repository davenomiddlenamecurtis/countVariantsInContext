#include "modelStrandBias.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
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

#define NCOMPONENT (2 + 2*4*2 + 2*16*3 + 2*64*2 + 2*256) // numbers of combinations for 1,2,3,4,5 bases

char chrNames[23][3], backgroundBaseTable[512][6], componentNames[4096][9], componentSequences[4096][7];
float backgroundCountTable[512][2];

int main(int argc, char* argv[])
{
	int i, c,f,l,comp,*useThis, *useThese,upToThree;
	float c1, c2,realBeta;
	FILE* fi,*fo;
	char *chr,line[1000],fn[100],seq[7], compPairBackgroundCountsFileRoot[100],component[6],intercept[20];
	glfModel model;
	useThis =(int*) calloc(NCOMPONENT,sizeof(int));
	useThese = (int*)calloc(NCOMPONENT, sizeof(int));
	for (comp = 0; comp < 2 + 2 * 4 * 2 + 2 * 16 * 3; ++comp)
		useThese[comp] = 1; // will include all combinations up to three in final models
	strcpy(intercept, "Intercept");
	// first, have a go at getting components
	l = 0;
	generateSequence(component, 1, 1);
	while (generateSequence(component, 1)) {
		if (strchr("CT", component[0])) {
			sprintf(componentSequences[l], "__%s__", component);
			sprintf(componentNames[l], "F1__%s__", component);
			++l;
		}
	}
	model.init(512, l);
	fi = fopen("compPairBackgroundCounts.total.txt", "r");
	i = 0;
	while (fgets(line, 999, fi)) {
		sscanf(line, "%s %*s %f %f", backgroundBaseTable[i], &c1, &c2);
		model.F[i] = c1 + c2;
		model.Y[i] = c1 / (c1 + c2);
		for (comp = 0; comp < l; ++comp)
			model.X[i][comp] = matches(backgroundBaseTable[i], componentSequences[comp], 5);
		++i;
	}
	for (comp = 0; comp < l; ++comp)
		model.name[comp] = componentNames[comp];
	model.name[comp] = intercept;
	fclose(fi);
	fo = fopen("modelBackgroundStrandBias.oneBase.txt", "w");
	evaluateModel(fo, &model, "backgroundStrandBiasOneBase");
	fclose(fo);
	generateSequence(component, 2, 1);
	while (generateSequence(component, 2)) {
		if (strchr("CT", component[1])) {
			sprintf(componentSequences[l], "_%s__", component);
			sprintf(componentNames[l], "U2_%s__", component);
			++l;
		}
		if (strchr("CT", component[0])) {
			sprintf(componentSequences[l], "__%s_", component);
			sprintf(componentNames[l], "D2__%s_", component);
			++l;
		}
	}
	for (c = 0; c < 22; ++c)
		sprintf(chrNames[c], "%d", c + 1);
	strcpy(chrNames[22], "X");
	strcpy(chrNames[23], "total");
	for (c = 0; c < 24; ++c) {
		model.init(512, l);
		chr = chrNames[c];
		sprintf(compPairBackgroundCountsFileRoot, "compPairBackgroundCounts.%s", chr);
		sprintf(fn, "%s.txt", compPairBackgroundCountsFileRoot);
		fi = fopen(fn, "r");
		i = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %*s %f %f", backgroundBaseTable[i], &c1, &c2);
			model.F[i] = c1 + c2;
			model.Y[i] = c1 / (c1 + c2);
			for (comp = 0; comp < l; ++comp)
				model.X[i][comp] = matches(backgroundBaseTable[i], componentSequences[comp],5);
			++i;
		}
		for (comp = 0; comp < l; ++comp)
			model.name[comp] = componentNames[comp];
		model.name[comp] = intercept;
		fclose(fi);
		sprintf(fn, "modelBackgroundStrandBias.%s.txt", chr);
		fo = fopen(fn, "w");
		evaluateModel(fo, &model, "backgroundStrandBias");
		fclose(fo);
	}
	generateSequence(component, 3, 1);
	while (generateSequence(component, 3)) {
		if (strchr("CT", component[2])) {
			sprintf(componentSequences[l], "%s__", component);
			sprintf(componentNames[l], "U3%s__", component);
			++l;
		}
		if (strchr("CT", component[1])) {
			sprintf(componentSequences[l], "_%s_", component);
			sprintf(componentNames[l], "F3_%s_", component);
			++l;
		}
		if (strchr("CT", component[0])) {
			sprintf(componentSequences[l], "__%s", component);
			sprintf(componentNames[l], "D3__%s", component);
			++l;
		}
	}
	model.init(512, l);
	fi = fopen("compPairBackgroundCounts.total.txt", "r");
	i = 0;
	while (fgets(line, 999, fi)) {
		sscanf(line, "%s %*s %f %f", backgroundBaseTable[i], &c1, &c2);
		model.F[i] = c1 + c2;
		model.Y[i] = c1 / (c1 + c2);
		for (comp = 0; comp < l; ++comp)
			model.X[i][comp] = matches(backgroundBaseTable[i], componentSequences[comp], 5);
		++i;
	}
	for (comp = 0; comp < l; ++comp)
		model.name[comp] = componentNames[comp];
	model.name[comp] = intercept;
	fclose(fi);
	fo = fopen("modelBackgroundStrandBias.trios.txt", "w");
	evaluateModel(fo, &model, "backgroundStrandBiasTrios");
	fclose(fo);

	upToThree = l;
	generateSequence(component, 4, 1);
	while (generateSequence(component, 4)) {
		if (strchr("CT", component[2])) {
			sprintf(componentSequences[l], "%s_", component);
			sprintf(componentNames[l], "U4%s_", component);
			++l;
		}
		if (strchr("CT", component[1])) {
			sprintf(componentSequences[l], "_%s", component);
			sprintf(componentNames[l], "D4_%s", component);
			++l;
		}
	}
	generateSequence(component, 5, 1);
	while (generateSequence(component, 5)) {
		if (strchr("CT", component[2])) {
			sprintf(componentSequences[l], "%s", component);
			sprintf(componentNames[l], "F5%s", component);
			++l;
		}
	}
	model.init(512, l);
	fi = fopen("compPairBackgroundCounts.total.txt", "r");
	i = 0;
	while (fgets(line, 999, fi)) {
		sscanf(line, "%s %*s %f %f", backgroundBaseTable[i], &c1, &c2);
		model.F[i] = c1 + c2;
		model.Y[i] = c1 / (c1 + c2);
		for (comp = 0; comp < l; ++comp)
			model.X[i][comp] = matches(backgroundBaseTable[i], componentSequences[comp], 5);
		++i;
	}
	for (comp = 0; comp < l; ++comp)
		model.name[comp] = componentNames[comp];
	model.name[comp] = intercept;
	fclose(fi);
	fo = fopen("modelBackgroundStrandBias.4and5.txt", "w");
	useThis[comp] = 1; // intercept
	for (comp = upToThree; comp < l; ++comp) {
		useThis[comp] = 1;
		evaluateModel(fo, &model, "backgroundStrandBias4and5", useThis);
		realBeta = model.SD[comp] ? model.beta[comp] / model.SD[comp] : 0.0;
		if (fabs(realBeta) >= 0.01)
			useThese[comp] = 1;
		useThis[comp] = 0;

	}
	fclose(fo);

	fo = fopen("modelBackgroundStrandBias.fullModel.txt", "w");
	useThese[comp] = 1; // intercept
	evaluateModel(fo, &model, "backgroundStrandBiasFullModel", useThese);
	fclose(fo);
	free(useThis);
	free(useThese);
	return 0;
}