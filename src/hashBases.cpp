#include <stdio.h>
#include "hashBases.hpp"

const char* baseNames = "CTGA";
const char* basePairs[2] = {
	"CTGA",
	"GACT"
};

baseHasher::baseHasher() {
	int i;
	for (i = 0; i < 4; ++i) {
		table[baseNames[i]] = i;
		bpTable[basePairs[0][i]] = basePairs[1][i];
	}
}

int baseHasher::getComplement(char* out, const char* in, int len) {
	int i;
	for (i = 0; i < len; ++i)
		out[len -1 - i] = bpTable[in[i]];
	return 1;
}

unsigned int baseHasher::hashBases(const char* in, int len)
{
	unsigned int rv=0;
	int i;
	for (i = 0; i < len; ++i)
	{
		rv |= table[*in++];
		rv <<= 2;
	}
	rv >>= 2;
	return rv;
}
// aim is to minimise array operations to maximise speed

int baseHasher::convertToSig(char* out, const char* in)
{
	sprintf(out, "%c[%c>%c]%c", in[1], in[2], in[5], in[3]);
	return 1;
}
