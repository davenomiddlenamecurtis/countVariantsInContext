#include "modelMutation.hpp"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if 0
in general, each row will have a sequence and a count 

4 REFs

4 ALTs

4x3 MUTs

For each MUT:

4 U2s
4 D2s

16 U3s
16 F3s
16 D3s

64 U4s
64 D4s

256 F5s

maybe first use F3s, then U3s, then D3s, then all 3s

#endif

int main(int argc, char* argv[])
{
	int i, c,f,l,b,bs,b1,b2,comp,*useThis, *useThese,upToThree;
	unsigned int hash;
	char target,source,first, second;
	float c1,realBeta;
	FILE* fi,*fo,*fplus,*fminus;
	char *chr,line[1000],fn[100],seq[7],component[6],intercept[20],minusStrandSeq[7],sig[8];
	baseHasher hasher;
	glfModel model;
	useThis =(int*) calloc(NCOMPONENTS,sizeof(int)); // numbers of combinations for 1,2,3,4,5 bases
	useThese = (int*)calloc(NCOMPONENTS, sizeof(int));
	fo = fopen("SBSVariants.txt", "w");
	generateSBSVariant(seq, 1);
	while (generateSBSVariant(seq)) {
		hasher.convertToSig(sig, seq);
		fprintf(fo, "%s\n", sig);
	}
	fclose(fo);

	strcpy(intercept, "Intercept");
	fi = fopen("backgroundCounts.total.txt", "r");
	i = 0;
	while (fgets(line, 999, fi)) {
		sscanf(line, "%s %f", backgroundSequenceTable[i], &c1);
		hash = hasher.hashBases(backgroundSequenceTable[i], 5);
		backgroundCountTable[hash] = c1;
		++i;
	}

	// for each base outcome, model background which leads to it
	sprintf(fn, "modelMutation.oneMut.txt");
	fo = fopen(fn, "w");
	for (b = 0; b < 4; ++b) {
		target = baseNames[b];
		l = 0;
		generateSequence(component, 1, 1);
		while (generateSequence(component, 1)) {
			if (component[0] != target) {
				sprintf(componentSequences[l], "__%c__%c", component[0], target);
				sprintf(componentNames[l], "M1__%c__%c", component[0], target);
				useThese[l] = 1;
				++l;
			}
		}
		model.init(NSEQUENCES / 4, l); // only those with correct target
		fi = fopen("counts.total.2.txt", "r");
		i = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
			if (sequenceTable[i][5] != target || sequenceTable[i][2] == target)
				continue;
			hash = hasher.hashBases(sequenceTable[i], 5); // just for background
			model.F[i] = backgroundCountTable[hash];
			model.Y[i] = c1 / backgroundCountTable[hash];
			for (comp = 0; comp < l; ++comp)
				model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
			++i;
		}
		for (comp = 0; comp < l; ++comp)
			model.name[comp] = componentNames[comp];
		model.name[comp] = intercept;
		fclose(fi);
		sprintf(line, "modelMutation.oneMut.%c", target);
		evaluateModel(fo, &model, line, useThese);
	}
	fclose(fo);

	getPredictedCounts("modelMutation.oneMut.txt","counts.total.2.txt","modelMutation.oneMut.predictions.txt");

	// for each base outcome, model background which leads to it
	sprintf(fn, "modelMutation.twoBases.txt");
	fo = fopen(fn, "w");
	for (b = 0; b < 4; ++b) {
		target = baseNames[b];
		l = 0;
		for (bs = 0; bs < 4; ++bs) {
			source = baseNames[bs];
			if (source == target)
				continue;
			generateSequence(component, 1, 1);
			while (generateSequence(component, 1)) {
				sprintf(componentSequences[l], "_%c%c__%c", component[0], source, target);
				sprintf(componentNames[l], "U2_%c%c__%c", component[0], source, target);
				useThese[l] = 1;
				++l;
				sprintf(componentSequences[l], "__%c%c_%c", source, component[0], target);
				sprintf(componentNames[l], "D2__%c%c_%c", source, component[0], target);
				useThese[l] = 1;
				++l;
			}
		}
		model.init(NSEQUENCES / 4, l); // only those with correct target
		fi = fopen("counts.total.2.txt", "r");
		i = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
			if (sequenceTable[i][5] != target || sequenceTable[i][2] == target)
				continue;
			hash = hasher.hashBases(sequenceTable[i], 5); // just for background
			model.F[i] = backgroundCountTable[hash];
			model.Y[i] = c1 / backgroundCountTable[hash];
			for (comp = 0; comp < l; ++comp)
				model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
			++i;
		}
		for (comp = 0; comp < l; ++comp)
			model.name[comp] = componentNames[comp];
		model.name[comp] = intercept;
		fclose(fi);
		sprintf(line, "modelMutation.twoBases.%c", target);
		evaluateModel(fo, &model, line, useThese);
	}
	fclose(fo);

	getPredictedCounts("modelMutation.twoBases.txt", "counts.total.2.txt", "modelMutation.twoBases.predictions.txt");

	sprintf(fn, "modelMutation.threeFlankingBases.txt");
	fplus = fopen("modelMutation.threeFlankingBases.plus.txt", "w");
	fminus = fopen("modelMutation.threeFlankingBases.minus.txt", "w");
	fclose(fplus);
	fclose(fminus);
	fo = fopen(fn, "w");
	for (b = 0; b < 4; ++b) {
		target = baseNames[b];
		l = 0;
		for (bs = 0; bs < 4; ++bs) {
			source = baseNames[bs];
			if (source == target)
				continue;
			for (b1 = 0; b1 < 4; ++b1) {
				first = baseNames[b1];
				for (b2 = 0; b2 < 4; ++b2) {
					second = baseNames[b2];
					sprintf(componentSequences[l], "_%c%c%c_%c", first, source, second, target);
					sprintf(componentNames[l], "F3_%c%c%c_%c", first, source, second, target);
					useThese[l] = 1;
					++l;
				}
			}
		}
		model.init(NSEQUENCES / 4, l); // only those with correct target
		fi = fopen("counts.total.2.txt", "r");
		i = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
			if (sequenceTable[i][5] != target || sequenceTable[i][2] == target)
				continue;
			hash = hasher.hashBases(sequenceTable[i], 5); // just for background
			model.F[i] = backgroundCountTable[hash];
			model.Y[i] = c1 / backgroundCountTable[hash];
			for (comp = 0; comp < l; ++comp)
				model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
			++i;
		}
		for (comp = 0; comp < l; ++comp)
			model.name[comp] = componentNames[comp];
		model.name[comp] = intercept;
		fclose(fi);
		sprintf(line, "modelMutation.threeFlankingBases.%c", target);
		evaluateModel(fo, &model, line, useThese);
		fplus = fopen("modelMutation.threeFlankingBases.plus.txt", "a");
		fminus = fopen("modelMutation.threeFlankingBases.minus.txt", "a");
		for (comp = 0; comp < l; ++comp) {
			strcpy(seq, model.name[comp]+2);
			realBeta= model.SD[comp] ? model.beta[comp] / model.SD[comp] : 0.0;
			if (seq[2]=='C'|| seq[2] == 'T') {
				hasher.convertToSig(sig, seq);
				fprintf(fplus, "%s\t%f\n", sig,realBeta);
			}
			else {
				hasher.getComplement(minusStrandSeq+1, seq+1, 3);
				hasher.getComplement(minusStrandSeq+5, seq+5, 1);
				minusStrandSeq[0] = '_';
				minusStrandSeq[4] = '_';
				minusStrandSeq[6] = '\0';
				hasher.convertToSig(sig, minusStrandSeq);
				fprintf(fminus, "%s\t%f\n", sig,realBeta);
			}
		}
		fclose(fplus);
		fclose(fminus);
	}
	fclose(fo);

	getPredictedCounts("modelMutation.threeFlankingBases.txt", "counts.total.2.txt", "modelMutation.threeFlankingBases.predictions.txt");

	sprintf(fn, "modelMutation.threeUpstreamBases.txt");
	fo = fopen(fn, "w");
	for (b = 0; b < 4; ++b) {
		target = baseNames[b];
		l = 0;
		for (bs = 0; bs < 4; ++bs) {
			source = baseNames[bs];
			if (source == target)
				continue;
			for (b1 = 0; b1 < 4; ++b1) {
				first = baseNames[b1];
				for (b2 = 0; b2 < 4; ++b2) {
					second = baseNames[b2];
					sprintf(componentSequences[l], "%c%c%c__%c", second, first, source, target);
					sprintf(componentNames[l], "U3%c%c%c__%c", second, first, source, target);
					useThese[l] = 1;
					++l;
				}
			}
		}
		model.init(NSEQUENCES / 4, l); // only those with correct target
		fi = fopen("counts.total.2.txt", "r");
		i = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
			if (sequenceTable[i][5] != target || sequenceTable[i][2] == target)
				continue;
			hash = hasher.hashBases(sequenceTable[i], 5); // just for background
			model.F[i] = backgroundCountTable[hash];
			model.Y[i] = c1 / backgroundCountTable[hash];
			for (comp = 0; comp < l; ++comp)
				model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
			++i;
		}
		for (comp = 0; comp < l; ++comp)
			model.name[comp] = componentNames[comp];
		model.name[comp] = intercept;
		fclose(fi);
		sprintf(line, "modelMutation.threeUpstreamBases.%c", target);
		evaluateModel(fo, &model, line, useThese);
	}
	fclose(fo);

	getPredictedCounts("modelMutation.threeUpstreamBases.txt", "counts.total.2.txt", "modelMutation.threeUpstreamBases.predictions.txt");


	sprintf(fn, "modelMutation.threeDownstreamBases.txt");
	fo = fopen(fn, "w");
	for (b = 0; b < 4; ++b) {
		target = baseNames[b];
		l = 0;
		for (bs = 0; bs < 4; ++bs) {
			source = baseNames[bs];
			if (source == target)
				continue;
			for (b1 = 0; b1 < 4; ++b1) {
				first = baseNames[b1];
				for (b2 = 0; b2 < 4; ++b2) {
					second = baseNames[b2];
					sprintf(componentSequences[l], "__%c%c%c%c", source, first, second, target);
					sprintf(componentNames[l], "D3__%c%c%c%c", source, first, second, target);
					useThese[l] = 1;
					++l;
				}
			}
	}
		model.init(NSEQUENCES / 4, l); // only those with correct target
		fi = fopen("counts.total.2.txt", "r");
		i = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
			if (sequenceTable[i][5] != target || sequenceTable[i][2] == target)
				continue;
			hash = hasher.hashBases(sequenceTable[i], 5); // just for background
			model.F[i] = backgroundCountTable[hash];
			model.Y[i] = c1 / backgroundCountTable[hash];
			for (comp = 0; comp < l; ++comp)
				model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
			++i;
		}
		for (comp = 0; comp < l; ++comp)
			model.name[comp] = componentNames[comp];
		model.name[comp] = intercept;
		fclose(fi);
		sprintf(line, "modelMutation.threeDownstreamBases.%c", target);
		evaluateModel(fo, &model, line, useThese);
	}
	fclose(fo);

	getPredictedCounts("modelMutation.threeDownstreamBases.txt", "counts.total.2.txt", "modelMutation.threeDownstreamBases.predictions.txt");


	sprintf(fn, "modelMutation.threeBases.txt");
	fo = fopen(fn, "w");
	for (b = 0; b < 4; ++b) {
		target = baseNames[b];
		l = 0;
		for (bs = 0; bs < 4; ++bs) {
			source = baseNames[bs];
			if (source == target)
				continue;
			for (b1 = 0; b1 < 4; ++b1) {
				first = baseNames[b1];
				for (b2 = 0; b2 < 4; ++b2) {
					second = baseNames[b2];
					sprintf(componentSequences[l], "_%c%c%c_%c", first, source, second, target);
					sprintf(componentNames[l], "F3_%c%c%c_%c", first, source, second, target);
					useThese[l] = 1;
					++l;
				}
			}
			for (b1 = 0; b1 < 4; ++b1) {
				first = baseNames[b1];
				for (b2 = 0; b2 < 4; ++b2) {
					second = baseNames[b2];
					sprintf(componentSequences[l], "%c%c%c__%c", second, first, source, target);
					sprintf(componentNames[l], "U3%c%c%c__%c", second, first, source, target);
					useThese[l] = 1;
					++l;
				}
			}
			for (b1 = 0; b1 < 4; ++b1) {
				first = baseNames[b1];
				for (b2 = 0; b2 < 4; ++b2) {
					second = baseNames[b2];
					sprintf(componentSequences[l], "__%c%c%c%c", source, first, second, target);
					sprintf(componentNames[l], "D3__%c%c%c%c", source, first, second, target);
					useThese[l] = 1;
					++l;
				}
			}
		}
		model.init(NSEQUENCES / 4, l); // only those with correct target
		fi = fopen("counts.total.2.txt", "r");
		i = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
			if (sequenceTable[i][5] != target || sequenceTable[i][2] == target)
				continue;
			hash = hasher.hashBases(sequenceTable[i], 5); // just for background
			model.F[i] = backgroundCountTable[hash];
			model.Y[i] = c1 / backgroundCountTable[hash];
			for (comp = 0; comp < l; ++comp)
				model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
			++i;
		}
		for (comp = 0; comp < l; ++comp)
			model.name[comp] = componentNames[comp];
		model.name[comp] = intercept;
		fclose(fi);
		sprintf(line, "modelMutation.threeBases.%c", target);
		evaluateModel(fo, &model, line, useThese);
	}
	fclose(fo);

	getPredictedCounts("modelMutation.threeBases.txt", "counts.total.2.txt", "modelMutation.threeBases.predictions.txt");


#if 0
	sprintf(fn, "modelMutation.twoAndThreeBases.txt");
	fo = fopen(fn, "w");
	for (b = 0; b < 4; ++b) {
		target = baseNames[b];
		l = 0;
		for (bs = 0; bs < 4; ++bs) {
			source = baseNames[bs];
			if (source == target)
				continue;
			generateSequence(component, 1, 1);
			while (generateSequence(component, 1)) {
				sprintf(componentSequences[l], "_%c%c__%c", component[0], source, target);
				sprintf(componentNames[l], "U2_%c%c__%c", component[0], source, target);
				useThese[l] = 1;
				++l;
				sprintf(componentSequences[l], "__%c%c_%c", source, component[0], target);
				sprintf(componentNames[l], "D2__%c%c_%c", source, component[0], target);
				useThese[l] = 1;
				++l;
			}
			generateSequence(component, 2, 1);
			while (generateSequence(component, 2)) {
				sprintf(componentSequences[l], "%c%c%c__%c", component[0], component[1], source, target);
				sprintf(componentNames[l], "U3%c%c%c__%c", component[0], component[1], source, target);
				useThese[l] = 1;
				++l;
				sprintf(componentSequences[l], "_%c%c%c_%c", component[0], source, component[1], target);
				sprintf(componentNames[l], "F3_%c%c%c_%c", component[0], source, component[1], target);
				useThese[l] = 1;
				++l;
				sprintf(componentSequences[l], "__%c%c%c%c", source, component[0], component[1], target);
				sprintf(componentNames[l], "D3__%c%c%c%c", source, component[0], component[1], target);
				useThese[l] = 1;
				++l;
			}
		}
		model.init(NSEQUENCES / 4, l); // only those with correct target
		fi = fopen("counts.total.2.txt", "r");
		i = 0;
		while (fgets(line, 999, fi)) {
			sscanf(line, "%s %f", sequenceTable[i], &c1); // no real need to keep this in a table, but leave it for now
			if (sequenceTable[i][5] != target || sequenceTable[i][2] == target)
				continue;
			hash = hasher.hashBases(sequenceTable[i], 5); // just for background
			model.F[i] = backgroundCountTable[hash];
			model.Y[i] = c1 / backgroundCountTable[hash];
			for (comp = 0; comp < l; ++comp)
				model.X[i][comp] = matches(sequenceTable[i], componentSequences[comp], 6);
			++i;
		}
		for (comp = 0; comp < l; ++comp)
			model.name[comp] = componentNames[comp];
		model.name[comp] = intercept;
		fclose(fi);
		sprintf(line, "modelMutation.twoAndThreeBases.%c", target);
		evaluateModel(fo, &model, line, useThese);
	}
	fclose(fo);

	getPredictedCounts("modelMutation.twoAndThreeBases.txt", "counts.total.2.txt", "modelMutation.twoAndThreeBases.predictions.txt");
#endif
	free(useThis);
	free(useThese);
	return 0;
}