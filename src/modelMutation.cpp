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
	int i, c,f,l,b,bs,b1,b2,comp,*useThis, *useThese,upToThree,numTrios;
	unsigned int hash;
	char target,source,first, second;
	float c1,realBeta,countTotal;
	double* keptBetas[4];
	FILE* fi,*fo,*fplus,*fminus,*flog,*fpairs,*fc;
	char *chr,line[1000],fn[100],seq[7],component[6],intercept[20],minusStrandSeq[7], sig[8],minusStrandSig[8];
	baseHasher hasher;
	glfModel model;
	useThis =(int*) calloc(NCOMPONENTS,sizeof(int)); // numbers of combinations for 1,2,3,4,5 bases
	useThese = (int*)calloc(NCOMPONENTS, sizeof(int));
	flog = fopen("modelMutation.log.txt", "w");
	fo = fopen("SBSVariantsForInputMatrix.txt", "w");
	generateSBSVariantForInputMatrix(seq, 1);
	while (generateSBSVariantForInputMatrix(seq)) {
		hasher.convertToSig(sig, seq);
		fprintf(fo, "%s\n", sig);
	}
	fclose(fo);
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
	countTotal = 0;
	while (fgets(line, 999, fi)) {
		sscanf(line, "%s %f", backgroundSequenceTable[i], &c1);
		hash = hasher.hashBases(backgroundSequenceTable[i], 5);
		backgroundCountTable[hash] = c1;
		countTotal += c1;
		++i;
	}
	meanBackgroundCount = countTotal / i;

	// for each base outcome, model background which leads to it
	sprintf(fn, "modelMutation.oneMut.txt");
	fo = fopen(fn, "w");
	sprintf(fn, "modelMutation.oneMut.fromCounts.txt");
	fc = fopen(fn, "w");
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
		// useThese[l] = 1; // do not use intercept, just ends up close to 0 anyway
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
		for (comp = 0; comp < l; ++comp)
			model.toFit[comp] = model.toUse[comp] = 1;
		setStartingBetasFromCounts(&model);
		printModel(fc, "oneMut from counts", 0, &model);
		sprintf(line, "modelMutation.oneMut.%c", target);
		evaluateModel(fo, &model, line, useThese);
	}
	fclose(fo);
	fclose(fc);

	getPredictedCounts(flog, "modelMutation.oneMut.txt", "correctedCounts.total.2.txt", "modelMutation.oneMut.predictions.txt");
	
	// for each base outcome, model background which leads to it
	sprintf(fn, "modelMutation.U2.txt");
	fo = fopen(fn, "w");
	sprintf(fn, "modelMutation.U2.fromCounts.txt");
	fc = fopen(fn, "w");
	fpairs = fopen("U2.complements.txt", "w"); // helpful for Excel to find complementary sequence
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
				hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
				fprintf(fpairs, "U2%s\tD2%s\n", componentSequences[l], minusStrandSeq);
				sprintf(componentNames[l], "U2_%c%c__%c", component[0], source, target);
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
		for (comp = 0; comp < l; ++comp)
			model.toFit[comp] = model.toUse[comp] = 1;
		setStartingBetasFromCounts(&model);
		printModel(fc, "U2 from counts", 0, &model);
		sprintf(line, "modelMutation.U2.%c", target);
		evaluateModel(fo, &model, line, useThese);
	}
	fclose(fo);
	fclose(fc);
	fclose(fpairs);

	getPredictedCounts(flog, "modelMutation.U2.txt", "correctedCounts.total.2.txt", "modelMutation.U2.predictions.txt");

	// for each base outcome, model background which leads to it
	sprintf(fn, "modelMutation.D2.txt");
	fo = fopen(fn, "w");
	sprintf(fn, "modelMutation.D2.fromCounts.txt");
	fc = fopen(fn, "w");
	fpairs = fopen("D2.complements.txt", "w");
	for (b = 0; b < 4; ++b) {
		target = baseNames[b];
		l = 0;
		for (bs = 0; bs < 4; ++bs) {
			source = baseNames[bs];
			if (source == target)
				continue;
			generateSequence(component, 1, 1);
			while (generateSequence(component, 1)) {
				sprintf(componentSequences[l], "__%c%c_%c", source, component[0], target);
				hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
				fprintf(fpairs, "D2%s\tU2%s\n", componentSequences[l], minusStrandSeq);
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
		for (comp = 0; comp < l; ++comp)
			model.toFit[comp] = model.toUse[comp] = 1;
		setStartingBetasFromCounts(&model);
		printModel(fc, "D2 from counts", 0, &model);
		sprintf(line, "modelMutation.D2.%c", target);
		evaluateModel(fo, &model, line, useThese);
	}
	fclose(fo);
	fclose(fc);

	getPredictedCounts(flog, "modelMutation.D2.txt", "correctedCounts.total.2.txt", "modelMutation.D2.predictions.txt");

	// for each base outcome, model background which leads to it
	sprintf(fn, "modelMutation.twoBases.txt");
	fo = fopen(fn, "w");
	sprintf(fn, "modelMutation.twoBases.fromCounts.txt");
	fc = fopen(fn, "w");
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
		for (comp = 0; comp < l; ++comp)
			model.toFit[comp] = model.toUse[comp] = 1;
		setStartingBetasFromCounts(&model);
		printModel(fc, "twoBases from counts", 0, &model);
		sprintf(line, "modelMutation.twoBases.%c", target);
		evaluateModel(fo, &model, line, useThese);
	}
	fclose(fo);
	fclose(fc);

	getPredictedCounts(flog, "modelMutation.twoBases.txt", "correctedCounts.total.2.txt", "modelMutation.twoBases.predictions.txt");

	sprintf(fn, "modelMutation.threeFlankingBases.txt");
	fplus = fopen("modelMutation.threeFlankingBases.plus.txt", "w");
	fminus = fopen("modelMutation.threeFlankingBases.minus.txt", "w");
	fclose(fplus);
	fclose(fminus);
	fo = fopen(fn, "w");
	sprintf(fn, "modelMutation.threeFlankingBases.fromCounts.txt");
	fc = fopen(fn, "w");
	fpairs = fopen("F3.complements.txt", "w"); // helpful for Excel to find complementary sequence
	for (b = 0; b < 4; ++b) {
		target = baseNames[b];
		keptBetas[b] = (double*)calloc(4 * 4 * 3 * 3+1, sizeof(double)); // enough for three sets of trios: F, U and D, plus 1 for the intercept
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
					hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
					fprintf(fpairs, "F3%s\tF3%s\n", componentSequences[l], minusStrandSeq);
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
		numTrios = comp; // will need to pull the intercept beta from this later
		fclose(fi);
		sprintf(line, "modelMutation.threeFlankingBases.%c", target);
		for (comp = 0; comp < l; ++comp)
			model.toFit[comp] = model.toUse[comp] = 1;
		setStartingBetasFromCounts(&model);
		printModel(fc, "F3 from counts",0, & model);
		evaluateModel(fo, &model, line, useThese);
		getBetas(keptBetas[b], &model);
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
#if 0
				hasher.getComplement(minusStrandSeq+1, seq+1, 3);
				hasher.getComplement(minusStrandSeq+5, seq+5, 1);
				minusStrandSeq[0] = '_';
				minusStrandSeq[4] = '_';
				minusStrandSeq[6] = '\0';
#else
				hasher.getVariantComplement(minusStrandSeq, seq);
#endif
				hasher.convertToSig(sig, minusStrandSeq);
				fprintf(fminus, "%s\t%f\n", sig,realBeta);
			}
		}
		fclose(fplus);
		fclose(fminus);
	}
	fclose(fo);
	fclose(fc);
	fclose(fpairs);

	getPredictedCounts(flog, "modelMutation.threeFlankingBases.txt", "correctedCounts.total.2.txt", "modelMutation.threeFlankingBases.predictions.txt");

	sprintf(fn, "modelMutation.threeUpstreamBases.txt");
	fo = fopen(fn, "w");
	sprintf(fn, "modelMutation.threeUpstreamBases.fromCounts.txt");
	fc = fopen(fn, "w");
	fpairs = fopen("U3.complements.txt", "w");
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
					hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
					fprintf(fpairs, "U3%s\tD3%s\n", componentSequences[l], minusStrandSeq);
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
		for (comp = 0; comp < l; ++comp)
			model.toFit[comp] = model.toUse[comp] = 1;
		setStartingBetasFromCounts(&model);
		printModel(fc, "threeUpstream from counts", 0, &model);
		sprintf(line, "modelMutation.threeUpstreamBases.%c", target);
		evaluateModel(fo, &model, line, useThese);
	}
	fclose(fo);
	fclose(fc);
	fclose(fpairs);

	getPredictedCounts(flog, "modelMutation.threeUpstreamBases.txt", "correctedCounts.total.2.txt", "modelMutation.threeUpstreamBases.predictions.txt");

	sprintf(fn, "modelMutation.threeDownstreamBases.txt");
	fo = fopen(fn, "w");
	sprintf(fn, "modelMutation.threeDownstreamBases.fromCounts.txt");
	fc = fopen(fn, "w");
	fpairs = fopen("D3.complements.txt", "w");
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
					hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
					fprintf(fpairs, "D3%s\tU3%s\n", componentSequences[l], minusStrandSeq);
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
		for (comp = 0; comp < l; ++comp)
			model.toFit[comp] = model.toUse[comp] = 1;
		setStartingBetasFromCounts(&model);
		printModel(fc, "threeDownstream from counts", 0, &model);
		sprintf(line, "modelMutation.threeDownstreamBases.%c", target);
		evaluateModel(fo, &model, line, useThese);
	}
	fclose(fo);
	fclose(fc);
	fclose(fpairs);

	getPredictedCounts(flog, "modelMutation.threeDownstreamBases.txt", "correctedCounts.total.2.txt", "modelMutation.threeDownstreamBases.predictions.txt");

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
				} // essential that these are generated in same order as previously so that beta match
			}
		}
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
		keptBetas[b][comp] = keptBetas[b][numTrios]; // deal with the intercept correctly, though at present not included
		// may change this to use betas from counts, which should be the same thing
		for (comp = numTrios; comp < l; ++comp)
			keptBetas[b][comp] = 0;
		sprintf(line, "modelMutation.threeBases.%c", target);
#if 1
		for (comp = 0; comp < l; ++comp)
			model.toFit[comp] = model.toUse[comp] = 1;
		setStartingBetasFromCounts(&model);
		evaluateModel(fo, &model, line, useThese, 0); // use starting betas from F3, U3 and D3 counts
#else
		evaluateModel(fo, &model, line, useThese,keptBetas[b]); // hope is that starting betas will achieve better maximisation
#endif
	}
	fclose(fo);

	getPredictedCounts(flog, "modelMutation.threeBases.txt", "correctedCounts.total.2.txt", "modelMutation.threeBases.predictions.txt");



	free(useThis);
	free(useThese);
	return 0;
}
