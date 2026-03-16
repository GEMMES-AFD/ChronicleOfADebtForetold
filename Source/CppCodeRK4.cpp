#include <iostream>
#include <fstream>
#include<math.h>
#include <Rcpp.h>
using namespace Rcpp; 

#define dim 73
#define dimIv 245
#define dimOut 391

//Note: All pieces of code beginning with a @ will be replaced by the required code by R before compiling
//For instance @AddDim will be replaced by the dimension of the model

// utility function that converts a Rcpp::List to a double**
// WARNING: do not forget to free the double** after use!
template<typename T>
void RcppListToPptr(Rcpp::List L, T**& pptr) {
	for (unsigned int it=0; it<L.size(); it++) {
		std::vector<double> tempVec = L[it];
		pptr[it] = (double*) malloc(sizeof(*pptr[it]) * tempVec.size());
		for (unsigned int it2=0; it2<tempVec.size(); it2++) {
			pptr[it][it2] = tempVec[it2];
		}
	}
}

void Func(double t, double* y, double* parms, double* ydot, double* x, double** dataExogVar, double** exogSamplingTime, int nExogVar, int* comptExogVar) {

ydot[71] = -parms[215] * y[71];
x[3] = parms[1] * y[0];
x[7] = y[31] + y[41];
x[12] = ((1.0/(1.0 + exp(parms[4] * t - parms[5]))) * (parms[6] - parms[6] * parms[7]) + parms[6] * parms[7]);
x[13] = (1.0/(1.0 + exp(-parms[209] * (t - parms[210])))) * parms[211];
x[22] = parms[221] + parms[222] * pow((y[62]/y[61]), parms[223]);
x[23] = (1.0/(1.0 + exp(parms[29] * t - parms[30]))) * (parms[11] - parms[11] * parms[12]) + parms[11] * parms[12];
x[24] = (1.0/(1.0 + exp(parms[29] * t - parms[30]))) * (parms[17] - parms[17] * parms[18]) + parms[17] * parms[18];
x[25] = (1.0/(1.0 + exp(parms[29] * t - parms[30]))) * (parms[23] - parms[23] * parms[24]) + parms[23] * parms[24];
x[26] = y[71] * y[65] * y[58] + y[9] * y[66] * y[64] * y[58];
x[29] = (1.0/(1.0 + exp(parms[38] * t - parms[39]))) * (parms[31] - parms[31] * parms[36] * parms[37]) + parms[31] * parms[36] * parms[37];
x[30] = parms[219] * pow((y[3]/y[2]), parms[220]);
x[34] = parms[40] - parms[41] * (y[1]/y[0] - parms[1]);
x[38] = (y[11] * (1.0 - y[7]) + y[63] * y[58] * (1.0 + parms[138]) * y[7]) * (1.0 + parms[142]);
x[40] = (y[11] * (1.0 - y[8]) + y[63] * y[58] * (1.0 + parms[212]) * y[8]) * (1.0 + parms[213]);
x[43] = y[63] * y[58]/y[11];
x[44] = y[63] * y[58];
x[47] = parms[224] * pow((y[3]/y[2]), parms[225]);
x[65] = parms[70] * y[35];
x[66] = parms[71] * y[17];
x[71] = ((1.0/(1.0 + exp(parms[74] * t - parms[75]))) * (parms[76] - parms[76] * parms[77]) + parms[76] * parms[77]);
x[78] = parms[84] * (y[17] + y[35] + y[15] * y[58]);
x[82] = 0.0;
x[83] = y[43] + y[34] + y[13];
x[97] = y[69] + parms[118];
x[98] = y[69] + parms[119];
x[100] = 0.155;
x[107] = 0.0;
x[121] = parms[145] * (y[71] * y[65] * y[58]);
x[129] = parms[150] * y[67];
x[131] = ((1.0/(1.0 + exp(parms[151] * t - parms[152]))) * (parms[153] - parms[153] * parms[154]) + parms[153] * parms[154]);
x[132] = parms[155] * y[42];
x[141] = y[30] + y[54];
x[151] = parms[179] * y[66] * y[63];
x[153] = (1.0/(1.0 + exp(parms[180] * t - parms[181]))) * (parms[182] - parms[182] * parms[183]) + parms[182] * parms[183];
x[163] = parms[193];
x[181] = parms[52] * y[59];
x[182] = parms[52] * y[61];
x[183] = parms[53] * y[67];
x[186] = parms[62] * (parms[61] * (y[15] + y[16]) - y[14]);
x[193] = parms[78] * y[18];
x[197] = parms[87] * (parms[86] * y[23] - y[22]);
x[214] = parms[128] * y[33];
x[216] = 0.0;
x[225] = parms[168] * (parms[167] * (y[49] + y[47]) - y[45]);
x[239] = parms[191] * y[60];
x[240] = parms[191] * y[62];
x[241] = parms[193] * y[63];
x[242] = parms[193] * y[64];
x[243] = parms[194] * y[65];
x[244] = parms[192] * y[66];
ydot[14] = x[186];
ydot[18] = x[193];
ydot[22] = x[197];
ydot[37] = x[214];
ydot[38] = x[216];
ydot[45] = x[225];
ydot[59] = x[181];
ydot[60] = x[239];
ydot[61] = x[182];
ydot[62] = x[240];
ydot[63] = x[241];
ydot[64] = x[242];
ydot[65] = x[243];
ydot[66] = x[244];
ydot[67] = x[183];
x[2] = parms[2] * (x[3] - y[1]);
x[14] = x[13] * parms[208] + y[72];
x[19] = x[23] * pow((y[11]/(y[63] * y[58] * (1.0 + parms[138]))), parms[13]) + parms[14] * pow((y[60]/y[59]), parms[15]);
x[20] = x[24] * pow((y[11]/(y[63] * y[58] * (1.0 + parms[138]))), parms[19]) + parms[20] * pow((y[60]/y[59]), parms[21]);
x[21] = x[25] * pow((y[11]/(y[63] * y[58] * (1.0 + parms[138]))), parms[25]) + parms[26] * pow((y[60]/y[59]), parms[27]);
x[28] = x[29] * pow((y[64] * y[58]/(y[11] * (1.0 + x[30]))), parms[32]) + parms[33] * (pow((y[59]/y[60]), parms[34]));
x[33] = (1.0 + x[34]) * y[10];
x[39] = x[26]/(y[71] + y[9] * y[66]);
x[64] = x[65] + x[66];
x[68] = parms[72] * y[32] * x[38];
x[69] = parms[73] * y[2] * x[38];
x[70] = x[71] * y[18];
x[79] = parms[85] * (x[78] - y[27]);
x[86] = parms[88] - parms[89]/(1.0 + exp(-parms[90] * (y[24]/x[83] - parms[91])));
x[99] = y[30] * y[24] + x[97] * y[57] * y[58] - x[100] * y[44];
x[130] = x[131] * x[129];
x[146] = (1.0 - parms[207]) * x[141];
x[152] = x[153] * y[4] * x[38];
x[177] = parms[214] * (x[22] - y[8]);
x[212] = y[33]/x[38] - parms[133] * y[32];
x[221] = y[40]/x[38] - parms[156] * y[42];
x[236] = x[186] + x[197] + x[225];
x[238] = parms[190] * (x[163] - y[69]);
ydot[8] = x[177];
ydot[32] = x[212];
ydot[42] = x[221];
ydot[56] = x[236];
ydot[69] = x[238];
x[1] = y[0] + x[2];
x[15] = x[14] * x[40];
x[27] = x[26]/x[39];
x[42] = parms[201] * (1.0 - parms[202]) * x[14] * x[40]/x[7];
x[67] = x[68] + x[69];
x[85] = y[30] - x[86];
x[150] = parms[200] * x[14] * x[40]/y[58];
x[154] = parms[184] * x[152];
x[157] = (1.0 - parms[184]) * x[152];
x[173] = x[14] - parms[10] * y[3];
x[174] = parms[16] * (x[19] - y[5]);
x[175] = parms[22] * (x[20] - y[6]);
x[176] = parms[28] * (x[21] - y[7]);
x[178] = parms[35] * (x[28] - y[9]);
x[180] = parms[43] * (x[33] - y[11]);
x[204] = x[79];
x[229] = parms[203] * parms[204] * x[14] * x[40]/y[58];
x[230] = parms[203] * (1.0 - parms[204]) * x[14] * x[40];
ydot[3] = x[173];
ydot[5] = x[174];
ydot[6] = x[175];
ydot[7] = x[176];
ydot[9] = x[178];
ydot[11] = x[180];
ydot[27] = x[204];
ydot[50] = x[229];
ydot[51] = x[230];
x[9] = parms[3] * x[1];
x[36] = (y[11] * (1.0 - y[5]) + y[63] * y[58] * (1.0 + parms[138]) * y[5]) * (1.0 + parms[139] + parms[140] + x[42]);
x[63] = x[67] + x[64];
x[88] = ((x[85] * y[43] + x[85] * (1.0 - parms[92]) * y[13] + x[85] * (1.0 - parms[93]) * y[34]) + y[30] * y[24])/(x[83] + y[24]);
x[96] = parms[114] + parms[115] * (x[180]/y[11] - parms[116]);
x[106] = (1.0/(1.0 + exp(-parms[123] * (x[85] - parms[124] - x[180]/y[11])))) * (parms[125] - parms[126]) + parms[126];
x[155] = parms[185] * x[154];
x[194] = (parms[81] * (x[181]/y[59]) + parms[82] * x[180]/y[11]) * y[20];
x[202] = (1.0 - parms[205]) * x[230];
x[222] = (parms[147] * (x[181]/y[59]) + parms[148] * x[180]/y[11]) * y[39];
x[234] = parms[205] * x[230];
ydot[20] = x[194];
ydot[39] = x[222];
ydot[52] = x[234];
ydot[53] = x[202];
x[0] = (y[4] + x[155]/x[38])/y[2] - parms[10];
x[41] = parms[201] * parms[202] * x[14] * x[40]/((y[11] * (1.0 - y[6]) + y[63] * y[58] * (1.0 + parms[138]) * y[6]) * (1.0 + parms[141]) * (x[9] + x[70] + x[130]));
x[72] = parms[79] * x[63];
x[87] = x[88] * (1.0 + y[28]);
x[117] = parms[139] * x[7] * ((1.0 - y[5]) * y[11] + y[5] * y[63] * y[58] * (1.0 + parms[138]))/x[36];
x[156] = x[154] - x[155];
x[172] = y[4] + x[155]/x[38] - parms[10] * y[2];
x[209] = parms[117] * (x[96] - y[30]);
ydot[2] = x[172];
ydot[30] = x[209];
x[10] = x[38] * y[4] + y[33] + x[72] + y[40] + x[155];
x[37] = (y[11] * (1.0 - y[6]) + y[63] * y[58] * (1.0 + parms[138]) * y[6]) * (1.0 + parms[141] + x[41]);
x[90] = x[87] * (1.0 + y[29]);
x[192] = x[72]/x[38] - parms[80] * y[19];
ydot[19] = x[192];
x[8] = x[37] * (x[9] + x[130] + x[70]);
x[75] = x[63] - x[37] * x[70] - parms[144] * x[63] - (1.0 + parms[83]) * y[20] * y[18];
x[127] = (1.0 + parms[149]) * y[39] * x[129] + x[37] * x[130] + parms[156] * x[38] * y[42];
x[5] = x[7] + x[8] + x[10] + x[15] + x[26];
x[6] = x[7]/x[36] + x[8]/x[37] + x[10]/x[38] + x[14] + x[26]/x[39];
x[18] = y[5] * (x[7]/x[36]) + y[6] * (x[8]/x[37]) + y[7] * (x[10]/x[38]) + y[8] * x[15]/x[38];
x[118] = parms[140] * x[7] * ((1.0 - y[5]) * y[11] + y[5] * y[63] * y[58] * (1.0 + parms[138]))/x[36] + parms[141] * x[8] * ((1.0 - y[6]) * y[11] + y[6] * y[63] * y[58] * (1.0 + parms[138]))/x[37] + parms[142] * x[10] * ((1.0 - y[7]) * y[11] + y[7] * y[63] * y[58] * (1.0 + parms[138]))/x[38];
x[119] = x[42] * x[7] * ((1.0 - y[5]) * y[11] + y[5] * y[63] * y[58] * (1.0 + parms[138]))/x[36] + x[41] * x[8] * ((1.0 - y[6]) * y[11] + y[6] * y[63] * y[58] * (1.0 + parms[138]))/x[37];
x[126] = x[127] + y[41];
x[4] = x[1] - x[18];
x[17] = x[18] * y[63] * y[58];
x[32] = x[7]/x[36] + x[10]/x[38] + x[14] + x[27] - x[18];
x[148] = x[26] - x[18] * y[63] * y[58];
x[169] = parms[0] * (x[6] - y[0]) + x[0] * y[0];
x[170] = x[1] - x[6];
x[208] = std::max(parms[120] * x[18] * y[63] - y[57], 0.0);
ydot[0] = x[169];
ydot[1] = x[170];
ydot[57] = x[208];
x[31] = y[31] + y[41] + x[65] + x[68] + x[127] + x[10] + x[15] + x[26] - x[17];
x[46] = x[4]/(y[59]);
x[89] = parms[94] + parms[95]/(1.0 + exp(-parms[96] * ((y[17] + y[15] * y[58] + y[16] * y[58])/y[11] * x[4])));
x[116] = parms[138] * x[17];
x[120] = parms[143] * y[11] * x[4] + parms[144] * x[63];
x[162] = parms[186] * pow((x[17]/(y[55] * y[58])), parms[187]);
x[164] = parms[195] * x[4] * y[11];
x[165] = parms[196] * x[4] * y[11];
x[166] = parms[197] * x[4] * y[11];
x[167] = parms[198] * x[4] * y[11];
x[168] = parms[199] * x[4] * y[11];
x[227] = -parms[172] * (y[70] * x[148]/y[58]);
x[228] = -(1.0 - parms[172]) * (y[70] * x[148]/y[58]);
x[233] = -parms[188] * x[148];
ydot[47] = x[227];
ydot[48] = x[233];
ydot[49] = x[228];
x[16] = (1.0/(1.0 + exp(-parms[216] * (t - parms[217])))) * parms[218] * std::max((0.013 * x[31] - x[15])/x[40], 0.0);
x[35] = ((1.0 + parms[51]) * y[12] * x[46] + x[37] * x[9] + parms[143] * y[11] * x[4])/x[4];
x[45] = x[46] + x[129] + y[18];
x[49] = x[5] - x[17] - x[116] - parms[143] * y[11] * x[4] - x[118] - x[117] - x[37] * x[9] - x[69] - x[66] - (1.0 + parms[51]) * y[12] * x[46];
x[92] = parms[103] + parms[104] * pow((x[162]), parms[102]);
x[102] = y[12] * x[46] + y[39] * x[129] + y[20] * y[18];
x[103] = parms[51] * y[12] * x[46] + parms[149] * y[39] * x[129] + parms[83] * y[20] * y[18];
x[128] = parms[158] * x[31];
x[133] = x[132] * x[38] + parms[226] * x[31];
x[134] = parms[160] * y[12] * (y[67] - x[129] - x[46] - y[18]) + parms[161] * y[12] * y[67];
x[142] = parms[173] + parms[174]/(exp(-parms[175] * ((y[46] + y[47] * y[58] + y[49] * y[58] + y[50] * y[58] + y[51])/x[31])));
x[160] = -(y[55] * y[58] + y[56] * y[58] - y[23] * y[58] - y[16] * y[58] - y[47] * y[58] - y[49] * y[58] - y[50] * y[58] - y[48] - y[52])/x[31];
x[161] = ((y[16] + y[23] + y[47] + y[49] + y[50] - y[56]) * y[58])/x[31];
x[185] = parms[60] * (parms[59] * y[12] * x[46] * (1.0 + parms[51]) - y[13]);
x[190] = 1.0/(1.0 + exp(-parms[65] * (x[162] - parms[68]))) * (parms[66] - parms[67]) + parms[67];
x[200] = 1.0/(1.0 + exp(-parms[109] * (x[162] - parms[112]))) * (parms[110] - parms[111]) + parms[111];
x[206] = parms[97] * (x[89] - y[28]);
ydot[13] = x[185];
ydot[28] = x[206];
ydot[72] = 2.0 * (x[16] - y[72]);
x[48] = 1.0 - x[45]/y[67];
x[50] = parms[44] * x[49];
x[51] = parms[45] * x[49];
x[52] = parms[46] * (x[49]);
x[73] = (1.0 - parms[162]) * x[134];
x[93] = y[69] + x[92];
x[95] = y[69] + parms[105] * x[92];
x[104] = x[103] + parms[121] * x[102];
x[135] = parms[162] * x[134];
x[143] = y[69] + parms[177] * x[92];
x[179] = parms[42] * (x[35] - y[10]);
x[184] = (parms[47] * (x[181]/y[59]) + parms[48] * (x[45]/y[67] - parms[49]) + parms[50] * x[180]/y[11]) * y[12];
x[217] = parms[134] * x[102];
x[219] = parms[159] * (x[128] - y[41]);
x[220] = parms[157] * (x[133] - y[40]);
x[232] = parms[176] * (x[142] - y[54]);
ydot[10] = x[179];
ydot[12] = x[184];
ydot[36] = x[217];
ydot[40] = x[220];
ydot[41] = x[219];
ydot[54] = x[232];
x[74] = (1.0 - parms[146]) * x[104];
x[94] = x[93] * (1.0 + parms[106] * y[28]);
x[110] = parms[129] - parms[130] * x[90] - parms[131] * x[48];
x[122] = parms[146] * x[104];
x[124] = x[126] + (1.0 + parms[149]) * y[39] * x[129] + x[37] * x[130] + y[40] + x[14] * x[40] + x[135];
x[144] = (1.0 - parms[178]) * x[143];
x[53] = x[49] - x[121] - x[87] * y[17] - x[94] * y[15] * y[58] - x[95] * y[16] * y[58] + x[85] * (1.0 - parms[92]) * y[13] - (x[52] + x[50] + x[51]);
x[76] = x[90] * y[35] + x[87] * y[17] + x[94] * y[15] * y[58] + x[141] * y[25] + y[53] * x[146] - (x[85] * y[43] + x[85] * (1.0 - parms[92]) * y[13] + x[85] * (1.0 - parms[93]) * y[34]) - x[93] * y[23] * y[58] + x[98] * y[21] * y[58] - y[30] * y[24] - x[37] * x[70] - parms[144] * x[63] - (1.0 + parms[83]) * y[20] * y[18] + x[74] - x[73] + x[67] + x[64];
x[145] = (1.0 - parms[206]) * x[144];
x[54] = (1.0 - parms[135]) * x[53];
x[77] = (1.0 - parms[136]) * x[76];
x[115] = parms[137] * x[102] + parms[135] * x[53] + parms[136] * x[76];
x[125] = x[141] * y[46] + x[143] * y[47] * y[58] + x[144] * y[49] * y[58] + y[50] * x[145] * y[58] + y[51] * x[146];
x[55] = x[54] - x[186] * y[58] - x[185];
x[56] = x[54]/(x[38] * y[2]);
x[80] = x[77] - x[79] - x[165];
x[114] = x[115] + x[116] + x[117] + x[118] + x[119] + x[120];
x[123] = x[124] + x[125];
x[11] = (x[12] + parms[8] * (x[56] - x[180]/y[11])) * y[2];
x[57] = (1.0 - parms[58]) * x[55];
x[61] = parms[58] * x[55] - x[164];
x[81] = x[80];
x[138] = parms[163] * x[123];
x[139] = parms[165] * x[123];
x[58] = std::max(0.0, (parms[54] + parms[55] * (y[71] * y[65])/(x[31]/y[58])) * x[57]);
x[59] = std::max(0.0, (parms[56] + parms[57] * (y[71] * y[65])/(x[31]/y[58])) * x[57]);
x[62] = x[38] * y[4] - x[61];
x[171] = parms[9] * (x[11] - y[4]);
x[223] = parms[164] * (x[138] - y[43]);
x[224] = parms[166] * (x[139] - y[44]);
ydot[4] = x[171];
ydot[43] = x[223];
ydot[44] = x[224];
x[60] = x[57] - x[58] - x[59];
x[113] = x[114] + x[121] + x[51] + x[122] + x[85] * y[43] + x[100] * y[44] + x[59] - x[166] + x[99];
x[149] = (x[151] + x[97] * y[57] + x[98] * y[21] + x[150] - x[145] * y[50] - x[143] * y[47] - x[144] * y[49] - x[93] * y[23] - x[95] * y[16]) * y[58] - x[58] - x[82] + x[168] - x[141] * y[48] - y[52] * x[146];
x[187] = parms[63] * (x[62]/y[58]);
x[189] = (1.0 - parms[69]) * parms[64] * (x[62]/y[58]);
ydot[16] = x[189];
x[101] = (1.0 - parms[137]) * x[102] + x[52] + x[103] + x[134] - x[104] - x[90] * y[35] + x[85] * (1.0 - parms[93]) * y[34] + x[60] + x[81] + x[151] * y[58] + x[50] - x[68] - x[65] + x[167];
x[136] = x[123] - x[113] - x[127];
x[147] = -(x[148] + x[149])/x[31];
x[188] = (1.0 - x[200]) * x[187];
x[198] = (parms[108] * x[204]/y[58] + x[187]);
ydot[15] = x[188];
x[91] = parms[98] + parms[99]/(1.0 + exp(-parms[100] * (y[35]/x[101])));
x[108] = y[38] * x[101];
x[109] = x[110] * x[101];
x[111] = x[101] - y[31];
x[137] = x[136] + x[223] + x[224] + x[225] * y[58];
x[140] = parms[169] + parms[170] * x[147];
x[191] = x[62] - x[188] * y[58] - x[189] * y[58] - x[156];
x[195] = x[198] - x[197] - x[188];
x[199] = (1.0 - x[200]) * x[198];
ydot[17] = x[191];
ydot[23] = x[199];
x[112] = y[33] - x[111];
x[158] = x[18] * y[63] + x[143] * y[47] + x[144] * y[49] + x[93] * y[23] + x[95] * y[16] + x[145] * y[50] + x[58]/y[58] + x[82]/y[58] + x[195] + x[236] + x[141] * y[48]/y[58] + y[52] * x[146]/y[58];
x[159] = x[26]/y[58] + x[151] + x[168]/y[58] + x[97] * y[57] + x[98] * y[21] + x[152]/y[58] + x[227] + x[228] + x[189] + x[199] + x[229] + x[150] + x[233]/y[58] + x[234]/y[58] - x[208];
x[207] = parms[101] * (x[91] - y[29]);
x[211] = parms[132] * (x[109] - y[33]);
x[215] = parms[127] * (x[108] - y[68]);
x[226] = parms[171] * (x[140] - y[70]);
x[231] = x[137] - x[227] * y[58] - x[228] * y[58] - x[229] * y[58] - x[230] - parms[200] * x[14] * x[40]/y[58];
x[235] = (x[148]/y[58] + x[149]/y[58] + x[152]/y[58] + x[227] + x[233]/y[58] + x[228] + x[189] + x[199] + x[229] + x[150] - x[236]);
ydot[29] = x[207];
ydot[33] = x[211];
ydot[46] = x[231];
ydot[55] = x[235];
ydot[68] = x[215];
ydot[70] = x[226];
x[105] = x[106] * x[101] + x[107] * (y[34] + y[36]) + x[215];
x[196] = x[235] - x[208];
x[201] = x[231] - x[233];
x[213] = x[215] + x[214];
x[237] = parms[189] * ((x[158] - x[159])/x[159]);
ydot[21] = x[196];
ydot[25] = x[201];
ydot[35] = x[213];
ydot[58] = x[237];
x[210] = parms[122] * (x[105] - y[31]);
x[218] = x[111] - y[33] + x[213] - x[217];
ydot[31] = x[210];
ydot[34] = x[218];
x[84] = (x[191] + x[213] + x[201] + x[202]) + x[72] + parms[113] * (x[223] + x[218] + x[185]) - (x[223] + x[218] + x[185] + x[204] + x[157]) - x[217] + (x[197] * y[58] + x[188] * y[58] - x[199] * y[58] + x[196] * y[58]);
x[203] = parms[113] * ((x[223] + x[218] + x[185]));
ydot[26] = x[203];
x[205] = std::max(x[84], -y[24]);
ydot[24] = x[205];
}
	
Rcpp::NumericMatrix RK4(int nt, 
                      double byT,
                      std::vector<double> Ry0,
                      std::vector<double> Rparms, 
                      double** dataExogVar,
                      double** exogSamplingTime, 
                      int nExogVar) {
	int it, it1;
	double *y = &Ry0[0];
	double *parms = &Rparms[0];
	double y1[dim], y2[dim], y3[dim], ydot0[dim], ydot1[dim], ydot2[dim], ydot3[dim], ydots[dim], x0[dimIv], x1[dimIv], x2[dimIv], x3[dimIv];
	Rcpp::NumericMatrix out(nt, dimOut);

	for (it=0; it<dim;it++) { //init out vector
		out(0, it)=y[it];
	}
	int comptExogVar[nExogVar];
	for (it=0; it<nExogVar; it++) comptExogVar[it]=1;

	// get intermediateVar and compute distance at t=0 //
	Func(0, y, parms, ydot0, x0, dataExogVar, exogSamplingTime, nExogVar, comptExogVar); 
						for (it1=0; it1<dim; it1++) {
							out(0, dim+it1) = ydot0[it1];
						}
						for (it1=0; it1<dimIv; it1++) {
							out(0, 2*dim+it1) = x0[it1];
						}
						 
	
	for (it=0; it<nExogVar; it++) comptExogVar[it]=1;
	
	for (it=0; it<(nt-1); it++) {

			if(it*byT>=4 && it*byT - 4<byT) { 
parms[215] = 0.025; 
} 

			

			Func(it*byT, y, parms, ydot0, x0, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);

			for (it1=0; it1<dim; it1++)
				y1[it1] = y[it1] + ydot0[it1]*0.5*byT;
			Func((it + 0.5)*byT, y1, parms, ydot1, x1, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);
			for (it1=0; it1<dim; it1++)
				y2[it1] = y[it1] + ydot1[it1]*0.5*byT;
			Func((it + 0.5)*byT, y2, parms, ydot2, x2, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);
			for (it1=0; it1<dim; it1++)
				y3[it1] = y[it1] + ydot2[it1]*byT;
			Func((it+1)*byT, y3, parms, ydot3, x3, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);
			for (it1=0; it1<dim; it1++) {
				ydots[it1] = (ydot0[it1] + 2.0*ydot1[it1] + 2.0*ydot2[it1] + ydot3[it1])/6.0;
			  out(it+1, it1) = y[it1];
				y[it1] = y[it1] + byT*ydots[it1];
				
			}
			
			for(it1=0;it1<dim;it1++){
							out(it+1, dim+it1) = ydots[it1];
						}
						for(it1=0;it1<dimIv;it1++){
							out(it+1, 2*dim+it1) = x0[it1];
						}
				
	}
	return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RK4(int nt, 
                        double byT,
                        std::vector<double> Ry0,
                        std::vector<double> Rparms, 
                        Rcpp::List RdataExogVar,
                        Rcpp::List RexogSamplingTime) {
	double** dataExogVar = (double**) malloc(sizeof(double*)*RdataExogVar.size());
	RcppListToPptr(RdataExogVar, dataExogVar);
	double** exogSamplingTime = (double**) malloc(sizeof(double*)*RexogSamplingTime.size());
	RcppListToPptr(RexogSamplingTime, exogSamplingTime);
	int nExogVar = RdataExogVar.size();
	Rcpp::NumericMatrix out = RK4(nt, byT, Ry0, Rparms, dataExogVar, exogSamplingTime, nExogVar);
	for (unsigned int it=0; it<RdataExogVar.size(); it++) {
		free(dataExogVar[it]);
		free(exogSamplingTime[it]);
	}
	free(dataExogVar);
	free(exogSamplingTime);
	
	return out;
}
