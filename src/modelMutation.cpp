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
	int i, c,f,l,b,bs,b1,b2,b3,b4,b5,comp,*useThis, *useThese,upToThree,numTrios;
	unsigned int hash;
	char target,source,first, second,third,fourth,fifth;
	float c1,realBeta,countTotal,LL,*fitted5LLs[4];
	double* keptBetas[4], *fitted5Betas[4],fitted5BetaThreshold=0.4;
	FILE* fi,*fo,*fplus,*fminus,*flog,*fpairs,*fc;
	char *chr,line[1000],fn[100],seq[7],component[6],intercept[20],minusStrandSeq[7], sig[8],minusStrandSig[8];
	baseHasher hasher;
	glfModel model;
	useThis =(int*) calloc(NCOMPONENTS,sizeof(int)); // numbers of combinations for 1,2,3,4,5 bases
	useThese = (int*)calloc(NCOMPONENTS, sizeof(int));
	// fGlfLog = fopen("modelMutationGfl.log.txt", "w");
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
		LL = model.getLnL();
		sprintf(line, "modelMutation.oneMut.%c", target);
		printModel(fo, line, LL, &model);
	}
	fclose(fo);

	getPredictedCounts(flog, "modelMutation.oneMut.txt", "correctedCounts.total.2.txt", "modelMutation.oneMut.predictions.txt");
	
	// for each base outcome, model background which leads to it
	sprintf(fn, "modelMutation.U2.txt");
	fo = fopen(fn, "w");
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
		LL = model.getLnL();
		sprintf(line, "modelMutation.U2.%c", target);
		printModel(fo, line, LL, &model);
	}
	fclose(fo);
	fclose(fpairs);

	getPredictedCounts(flog, "modelMutation.U2.txt", "correctedCounts.total.2.txt", "modelMutation.U2.predictions.txt");

	// for each base outcome, model background which leads to it
	sprintf(fn, "modelMutation.D2.txt");
	fo = fopen(fn, "w");
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
		LL = model.getLnL();
		sprintf(line, "modelMutation.D2.%c", target);
		printModel(fo, line, LL, &model);
	}
	fclose(fo);
	fclose(fpairs);

	getPredictedCounts(flog, "modelMutation.D2.txt", "correctedCounts.total.2.txt", "modelMutation.D2.predictions.txt");

	sprintf(fn, "modelMutation.threeFlankingBases.txt");
	fplus = fopen("modelMutation.threeFlankingBases.plus.txt", "w");
	fminus = fopen("modelMutation.threeFlankingBases.minus.txt", "w");
	fclose(fplus);
	fclose(fminus);
	fo = fopen(fn, "w");
	fpairs = fopen("F3.complements.txt", "w"); // helpful for Excel to find complementary sequence
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
		LL = model.getLnL();
		printModel(fo, line,LL, & model);
		fplus = fopen("modelMutation.threeFlankingBases.plus.txt", "a");
		fminus = fopen("modelMutation.threeFlankingBases.minus.txt", "a");
		for (comp = 0; comp < l; ++comp) {
			strcpy(seq, model.name[comp]+2);
			realBeta= model.beta[comp];
			if (seq[2]=='C'|| seq[2] == 'T') {
				hasher.convertToSig(sig, seq);
				fprintf(fplus, "%s\t%f\n", sig,realBeta);
			}
			else {
				hasher.getVariantComplement(minusStrandSeq, seq);
				hasher.convertToSig(sig, minusStrandSeq);
				fprintf(fminus, "%s\t%f\n", sig,realBeta);
			}
		}
		fclose(fplus);
		fclose(fminus);
	}
	fclose(fo);
	fclose(fpairs);

	getPredictedCounts(flog, "modelMutation.threeFlankingBases.txt", "correctedCounts.total.2.txt", "modelMutation.threeFlankingBases.predictions.txt");

	sprintf(fn, "modelMutation.threeUpstreamBases.txt");
	fo = fopen(fn, "w");
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
		LL = model.getLnL();
		sprintf(line, "modelMutation.threeUpstreamBases.%c", target);
		printModel(fo, line, LL, &model);
	}
	fclose(fo);
	fclose(fpairs);

	getPredictedCounts(flog, "modelMutation.threeUpstreamBases.txt", "correctedCounts.total.2.txt", "modelMutation.threeUpstreamBases.predictions.txt");

	sprintf(fn, "modelMutation.threeDownstreamBases.txt");
	fo = fopen(fn, "w");
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
		LL = model.getLnL();
		sprintf(line, "modelMutation.threeDownstreamBases.%c", target);
		printModel(fo, line, LL, &model);
	}
	fclose(fo);
	fclose(fpairs);

	getPredictedCounts(flog, "modelMutation.threeDownstreamBases.txt", "correctedCounts.total.2.txt", "modelMutation.threeDownstreamBases.predictions.txt");

	sprintf(fn, "modelMutation.fourUpstreamBases.txt");
	fo = fopen(fn, "w");
	fpairs = fopen("U4.complements.txt", "w");
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
					for (b3 = 0; b3 < 4; ++b3) {
						third = baseNames[b3];
						sprintf(componentSequences[l], "%c%c%c%c_%c", third, second, source, first, target);
						hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
						fprintf(fpairs, "U4%s\tD4%s\n", componentSequences[l], minusStrandSeq);
						sprintf(componentNames[l], "U4%s", componentSequences[l]);
						useThese[l] = 1;
						++l;
					}
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
		LL = model.getLnL();
		sprintf(line, "modelMutation.fourUpstreamBases.%c", target);
		printModel(fo, line, LL, &model);
	}
	fclose(fo);
	fclose(fpairs);

	getPredictedCounts(flog, "modelMutation.fourUpstreamBases.txt", "correctedCounts.total.2.txt", "modelMutation.fourUpstreamBases.predictions.txt");

	sprintf(fn, "modelMutation.fourDownstreamBases.txt");
	fo = fopen(fn, "w");
	fpairs = fopen("D4.complements.txt", "w");
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
					for (b3 = 0; b3 < 4; ++b3) {
						third = baseNames[b3];
						sprintf(componentSequences[l], "_%c%c%c%c%c", first, source, second, third, target);
						hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
						fprintf(fpairs, "D3%s\tU3%s\n", componentSequences[l], minusStrandSeq);
						sprintf(componentNames[l], "D4%s", componentSequences[l]);
						useThese[l] = 1;
						++l;
					}
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
		LL = model.getLnL();
		sprintf(line, "modelMutation.fourDownstreamBases.%c", target);
		printModel(fo, line, LL, &model);
	}
	fclose(fo);
	fclose(fpairs);

	getPredictedCounts(flog, "modelMutation.fourDownstreamBases.txt", "correctedCounts.total.2.txt", "modelMutation.fourDownstreamBases.predictions.txt");

	sprintf(fn, "modelMutation.five.txt");
	fo = fopen(fn, "w");
	fclose(fo);
	sprintf(fn, "modelMutation.fullModel.txt");
	fo = fopen(fn, "w");
	fclose(fo);
	fo=fopen("fitted5Betas.txt", "w");
	fclose(fo);
	fpairs = fopen("F5.complements.txt", "w");
	for (b = 0; b < 4; ++b) {
		fitted5Betas[b] = (double*)calloc(NCOMPONENTS, sizeof(double));
		fitted5LLs[b] = (float*)calloc(NCOMPONENTS, sizeof(float));
	}
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
		upToThree = l;

			for (bs = 0; bs < 4; ++bs) {
				source = baseNames[bs];
				if (source == target)
					continue;
				for (b1 = 0; b1 < 4; ++b1) {
					first = baseNames[b1];
					for (b2 = 0; b2 < 4; ++b2) {
						second = baseNames[b2];
						for (b3 = 0; b3 < 4; ++b3) {
							third = baseNames[b3];
							for (b4 = 0; b4 < 4; ++b4) {
								fourth = baseNames[b4];
								sprintf(componentSequences[l], "%c%c%c%c%c%c", first, second, source, third, fourth, target);
								hasher.getVariantComplement(minusStrandSeq, componentSequences[l]);
								fprintf(fpairs, "F5%s\tF5%s\n", componentSequences[l], minusStrandSeq);
								sprintf(componentNames[l], "F5%s", componentSequences[l]);
								++l;
							}
						}
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

		for (comp=0;comp<upToThree;++comp)
			useThis[comp] = 1;
		for (comp = upToThree; comp < l; ++comp)
			useThis[comp] = 0;
		for (c = 0; c < l; ++c)
			model.toFit[c] = model.toUse[c] = useThis[c];
		fo = fopen("modelMutation.F5.txt", "a");
		for (comp = upToThree; comp < l; ++comp) {
			setStartingBetasFromCounts(&model);
			model.toFit[comp] = model.toUse[comp] = useThis[comp]=1;
			evaluateModel(fo, &model, "modelMutation.F5",useThis);
			realBeta = model.beta[comp];
			fitted5Betas[b][comp] = realBeta;
			fitted5LLs[b][comp] = model.getLnL();
			model.toFit[comp] = model.toUse[comp] = useThis[comp] = 0;
		}
		fclose(fo);
		fo=fopen("fitted5Betas.txt", "a");
		for (comp = upToThree; comp < l; ++comp)
			fprintf(fo, "%s\t%f\t%f\n", componentNames[comp],fitted5Betas[b][comp],fitted5LLs[b][comp]);
		fclose(fo);

		sprintf(fn, "modelMutation.fullModel.txt");
		fo = fopen(fn,"a");
		for (comp = upToThree; comp < l; ++comp)
			if (fabs(fitted5Betas[b][comp]) > fitted5BetaThreshold)
				useThese[comp] = 1;
		for (comp = 0; comp < upToThree; ++comp)
			model.toFit[comp] = model.toUse[comp] = 1;
		setStartingBetasFromCounts(&model);

		sprintf(line, "modelMutation.fullModel.%c", target);
		for (comp = upToThree; comp < l; ++comp)
			model.toFit[comp] = model.toUse[comp] = useThese[l];
		evaluateModel(fo, &model, line, useThese, 0); // use starting betas from F3counts
		fclose(fo);
	}
	fclose(fpairs);

	getPredictedCounts(flog, "modelMutation.fullModel.txt", "correctedCounts.total.2.txt", "modelMutation.fullModel.predictions.txt");

	free(useThis);
	free(useThese);
	// fclose(fGlfLog);
	return 0;
}
