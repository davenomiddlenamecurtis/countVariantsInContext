#ifndef RUNMODELSHPP
#define RUNMODELSHPP
#include <stdio.h>
#include "dcerror.hpp"
#include "glfModel.hpp"

int generateSequence(char* s, int len, int start = 0);
int generateSBSVariant(char* s, int start=0);
int generateSBSVariantForInputMatrix(char* s, int start=0);
int matches(char* seq, char* comp, int len);
void printModel(FILE* fo, const char* LLstr, double LL, glfModel* m);
float evaluateModel(FILE* fo, glfModel* m, const char* name, int* useThese = 0,double *startingBetas=0,int useLinearRegression=0);
int getBetas(double* beta, glfModel* m);
void setStartingBetasFromCounts(glfModel *m);

#endif

