#include <iostream>
#include <fstream>
#include<math.h>
#include <Rcpp.h>
using namespace Rcpp; 

#define dim 78
#define dimIv 271
#define dimOut 427

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

ydot[72] = -parms[220] * y[72];
x[3] = parms[1] * y[0];
x[7] = y[31] + y[41];
x[12] = ((1.0/(1.0 + exp(parms[4] * t - parms[5]))) * (parms[6] - parms[6] * parms[7]) + parms[6] * parms[7]);
x[13] = parms[213] * (parms[227] + parms[224] * pow(std::max(0.0000001, (t - 4.0)), parms[226]) * exp(-parms[225] * (t - 4.0)));
x[21] = parms[245] + parms[246] * pow((y[62]/y[61]), parms[247]);
x[23] = (1.0/(1.0 + exp(parms[29] * t - parms[30]))) * (parms[11] - parms[11] * parms[12]) + parms[11] * parms[12];
x[24] = (1.0/(1.0 + exp(parms[29] * t - parms[30]))) * (parms[17] - parms[17] * parms[18]) + parms[17] * parms[18];
x[25] = (1.0/(1.0 + exp(parms[29] * t - parms[30]))) * (parms[23] - parms[23] * parms[24]) + parms[23] * parms[24];
x[26] = y[72] * y[65] * y[58] + y[9] * y[66] * y[64] * y[58];
x[29] = (1.0/(1.0 + exp(parms[37] * t - parms[38]))) * (parms[31] - parms[31] * parms[36]) + parms[31] * parms[36];
x[30] = parms[243] * pow((y[3]/y[2]), parms[244]);
x[34] = parms[39] - parms[40] * (y[1]/y[0] - parms[1]);
x[38] = (y[11] * (1.0 - y[7]) + y[63] * y[58] * (1.0 + y[74]) * y[7]) * (1.0 + parms[142]);
x[40] = (y[11] * (1.0 - y[8]) + y[63] * y[58] * (1.0 + parms[217]) * y[8]) * (1.0 + parms[218]);
x[43] = y[63] * y[58]/y[11];
x[44] = y[63] * y[58];
x[47] = parms[248] * pow((y[3]/y[2]), parms[249]);
x[65] = parms[68] * y[35];
x[66] = parms[69] * y[17];
x[71] = ((1.0/(1.0 + exp(parms[72] * t - parms[73]))) * (parms[74] - parms[74] * parms[75]) + parms[74] * parms[75]);
x[78] = parms[82] * (y[17] + y[35] + y[15] * y[58]);
x[82] = 0.0;
x[83] = y[43] + y[34] + y[13];
x[97] = y[70] + parms[115];
x[98] = y[70] + parms[116];
x[100] = 0.155;
x[107] = 0.0;
x[124] = parms[145] * (y[72] * y[65] * y[58]);
x[135] = parms[150] * y[67];
x[137] = ((1.0/(1.0 + exp(parms[151] * t - parms[152]))) * (parms[153] - parms[153] * parms[154]) + parms[153] * parms[154]);
x[150] = y[30] + y[54];
x[160] = parms[184] * y[66] * y[63];
x[162] = (1.0/(1.0 + exp(parms[185] * t - parms[186]))) * (parms[187] - parms[187] * parms[188]) + parms[187] * parms[188];
x[172] = parms[198];
x[203] = parms[51] * y[59];
x[204] = parms[51] * y[61];
x[205] = parms[52] * y[67];
x[208] = parms[61] * (parms[60] * (y[15] + y[16]) - y[14]);
x[215] = parms[76] * y[18];
x[219] = parms[85] * (parms[84] * y[23] - y[22]);
x[236] = parms[125] * y[33];
x[238] = 0.0;
x[247] = parms[173] * (parms[172] * (y[49] + y[75] + y[47] + y[77]) - y[45]);
x[264] = parms[196] * y[60];
x[265] = parms[196] * y[62];
x[266] = parms[198] * y[63];
x[267] = parms[198] * y[64];
x[268] = parms[199] * y[65];
x[269] = parms[197] * y[66];
x[270] = y[75];
ydot[14] = x[208];
ydot[18] = x[215];
ydot[22] = x[219];
ydot[37] = x[236];
ydot[38] = x[238];
ydot[45] = x[247];
ydot[59] = x[203];
ydot[60] = x[264];
ydot[61] = x[204];
ydot[62] = x[265];
ydot[63] = x[266];
ydot[64] = x[267];
ydot[65] = x[268];
ydot[66] = x[269];
ydot[67] = x[205];
x[2] = parms[2] * (x[3] - y[1]);
x[14] = x[13] * x[40];
x[18] = x[23] * pow((y[11]/(y[63] * y[58] * (1.0 + y[74]))), parms[13]) + parms[14] * pow((y[60]/y[59]), parms[15]);
x[19] = x[24] * pow((y[11]/(y[63] * y[58] * (1.0 + y[74]))), parms[19]) + parms[20] * pow((y[60]/y[59]), parms[21]);
x[20] = x[25] * pow((y[11]/(y[63] * y[58] * (1.0 + y[74]))), parms[25]) + parms[26] * pow((y[60]/y[59]), parms[27]);
x[28] = x[29] * pow((y[64] * y[58]/(y[11] * (1.0 + x[30]))), parms[32]) + parms[33] * (pow((y[59]/y[60]), parms[34]));
x[33] = (1.0 + x[34]) * y[10];
x[39] = x[26]/(y[72] + y[9] * y[66]);
x[42] = parms[206] * (1.0 - parms[207]) * x[13] * x[40]/x[7];
x[64] = x[65] + x[66];
x[68] = parms[70] * y[32] * x[38];
x[69] = parms[71] * y[2] * x[38];
x[70] = x[71] * y[18];
x[79] = parms[83] * (x[78] - y[27]);
x[86] = parms[86] - parms[87]/(1.0 + exp(-parms[88] * (y[24]/x[83] - parms[89])));
x[99] = y[30] * y[24] + x[97] * y[57] * y[58] - x[100] * y[44];
x[136] = x[137] * x[135];
x[155] = (1.0 - parms[212]) * x[150];
x[159] = parms[205] * x[13] * x[40]/y[58];
x[161] = x[162] * y[4] * x[38];
x[195] = x[13] - parms[10] * y[3];
x[199] = parms[219] * (x[21] - y[8]);
x[234] = y[33]/x[38] - parms[130] * y[32];
x[243] = y[40]/x[38] - parms[161] * y[42];
x[252] = parms[208] * (1.0 - parms[209]) * x[13] * x[40];
x[261] = x[208] + x[219] + x[247];
x[263] = parms[195] * (x[172] - y[70]);
ydot[3] = x[195];
ydot[8] = x[199];
ydot[32] = x[234];
ydot[42] = x[243];
ydot[51] = x[252];
ydot[56] = x[261];
ydot[70] = x[263];
x[1] = y[0] + x[2];
x[27] = x[26]/x[39];
x[36] = (y[11] * (1.0 - y[5]) + y[63] * y[58] * (1.0 + y[74]) * y[5]) * (1.0 + parms[139] + parms[140] + x[42]);
x[67] = x[68] + x[69];
x[85] = y[30] - x[86];
x[163] = parms[189] * x[161];
x[166] = (1.0 - parms[189]) * x[161];
x[196] = parms[16] * (x[18] - y[5]);
x[197] = parms[22] * (x[19] - y[6]);
x[198] = parms[28] * (x[20] - y[7]);
x[200] = parms[35] * (x[28] - y[9]);
x[202] = parms[42] * (x[33] - y[11]);
x[224] = (1.0 - parms[210]) * x[252];
x[226] = x[79];
x[259] = parms[210] * x[252];
ydot[5] = x[196];
ydot[6] = x[197];
ydot[7] = x[198];
ydot[9] = x[200];
ydot[11] = x[202];
ydot[27] = x[226];
ydot[52] = x[259];
ydot[53] = x[224];
x[9] = parms[3] * x[1];
x[63] = x[67] + x[64];
x[88] = ((x[85] * y[43] + x[85] * (1.0 - parms[90]) * y[13] + x[85] * (1.0 - parms[91]) * y[34]) + y[30] * y[24])/(x[83] + y[24]);
x[96] = parms[111] + parms[112] * (x[202]/y[11] - parms[113]);
x[106] = (1.0/(1.0 + exp(-parms[120] * (x[85] - parms[121] - x[202]/y[11])))) * (parms[122] - parms[123]) + parms[123];
x[120] = parms[139] * x[7] * ((1.0 - y[5]) * y[11] + y[5] * y[63] * y[58] * (1.0 + y[74]))/x[36];
x[164] = parms[190] * x[163];
x[216] = (parms[79] * (x[203]/y[59]) + parms[80] * x[202]/y[11]) * y[20];
ydot[20] = x[216];
x[0] = (y[4] + x[164]/x[38])/y[2] - parms[10];
x[41] = parms[206] * parms[207] * x[13] * x[40]/((y[11] * (1.0 - y[6]) + y[63] * y[58] * (1.0 + y[74]) * y[6]) * (1.0 + parms[141]) * (x[9] + x[70] + x[136]));
x[72] = parms[77] * x[63];
x[87] = x[88] * (1.0 + y[28]);
x[165] = x[163] - x[164];
x[174] = ydot[11]/y[11];
x[194] = y[4] + x[164]/x[38] - parms[10] * y[2];
x[231] = parms[114] * (x[96] - y[30]);
ydot[2] = x[194];
ydot[30] = x[231];
x[10] = x[38] * y[4] + y[33] + x[72] + y[40] + x[164];
x[37] = (y[11] * (1.0 - y[6]) + y[63] * y[58] * (1.0 + y[74]) * y[6]) * (1.0 + parms[141] + x[41]);
x[90] = x[87] * (1.0 + y[29]);
x[214] = x[72]/x[38] - parms[78] * y[19];
ydot[19] = x[214];
x[8] = x[37] * (x[9] + x[136] + x[70]);
x[75] = x[63] - x[37] * x[70] - parms[144] * x[63] - (1.0 + parms[81]) * y[20] * y[18];
x[133] = (1.0 + parms[149]) * y[39] * x[135] + x[37] * x[136] + parms[161] * x[38] * y[42];
x[5] = x[7] + x[8] + x[10] + x[14] + x[26];
x[6] = x[7]/x[36] + x[8]/x[37] + x[10]/x[38] + x[13] + x[26]/x[39];
x[17] = y[5] * (x[7]/x[36]) + y[6] * (x[8]/x[37]) + y[7] * (x[10]/x[38]) + y[8] * x[14]/x[38];
x[121] = parms[140] * x[7] * ((1.0 - y[5]) * y[11] + y[5] * y[63] * y[58] * (1.0 + y[74]))/x[36] + parms[141] * x[8] * ((1.0 - y[6]) * y[11] + y[6] * y[63] * y[58] * (1.0 + y[74]))/x[37] + parms[142] * x[10] * ((1.0 - y[7]) * y[11] + y[7] * y[63] * y[58] * (1.0 + y[74]))/x[38];
x[122] = x[42] * x[7] * ((1.0 - y[5]) * y[11] + y[5] * y[63] * y[58] * (1.0 + y[74]))/x[36] + x[41] * x[8] * ((1.0 - y[6]) * y[11] + y[6] * y[63] * y[58] * (1.0 + y[74]))/x[37];
x[132] = x[133] + y[41];
x[4] = x[1] - x[17];
x[16] = x[17] * y[63] * y[58];
x[32] = x[7]/x[36] + x[10]/x[38] + x[13] + x[27] - x[17];
x[157] = x[26] - x[17] * y[63] * y[58];
x[191] = parms[0] * (x[6] - y[0]) + x[0] * y[0];
x[192] = x[1] - x[6];
ydot[0] = x[191];
ydot[1] = x[192];
x[31] = y[31] + y[41] + x[65] + x[68] + x[133] + x[10] + x[14] + x[26] - x[16];
x[46] = x[4]/(y[59]);
x[89] = parms[92] + parms[93]/(1.0 + exp(-parms[94] * ((y[17] + y[15] * y[58] + y[16] * y[58])/y[11] * x[4])));
x[119] = y[74] * x[16];
x[123] = parms[143] * y[11] * x[4] + parms[144] * x[63];
x[171] = parms[191] * pow((x[16]/(y[55] * y[58])), parms[192]);
x[182] = parms[200] * x[4] * y[11];
x[183] = parms[201] * x[4] * y[11];
x[184] = parms[202] * x[4] * y[11];
x[185] = parms[203] * x[4] * y[11];
x[186] = parms[204] * x[4] * y[11];
x[258] = -parms[193] * x[157];
ydot[48] = x[258];
x[15] = (1.0/(1.0 + exp(-parms[221] * (t - parms[222])))) * parms[223] * std::max((0.013 * x[31] - x[14])/x[40], 0.0);
x[35] = ((1.0 + parms[50]) * y[12] * x[46] + x[37] * x[9] + parms[143] * y[11] * x[4])/x[4];
x[45] = x[46] + x[135] + y[18];
x[49] = x[5] - x[16] - x[119] - parms[143] * y[11] * x[4] - x[121] - x[120] - x[37] * x[9] - x[69] - x[66] - (1.0 + parms[50]) * y[12] * x[46];
x[92] = (parms[101] + parms[102] * pow((x[171]), parms[100])) * (1.0 + parms[261] * parms[260] * parms[257] * (t >= 4.0) * exp(-sqrt(parms[256]) * pow((t - 0.5 * (4.0 + 2.0/parms[256])), 2.0)));
x[102] = y[12] * x[46] + y[39] * x[135] + y[20] * y[18];
x[103] = parms[50] * y[12] * x[46] + parms[149] * y[39] * x[135] + parms[81] * y[20] * y[18];
x[139] = 1.0/(1.0 + exp(-parms[158] * ((y[46] + (y[47] + y[77] * (1.0 - parms[252])) * y[58] + (y[49] + y[75] * (1.0 - parms[252])) * y[58] + (y[50] + y[76] * (1.0 - parms[252])) * y[58] + y[51])/x[31]) - parms[160]));
x[142] = (parms[165] * y[12] * (y[68] - x[135] - x[46] - y[18]) + parms[166] * y[12] * y[67]);
x[145] = parms[159] * x[31];
x[151] = parms[178] + parms[179]/(exp(-parms[180] * ((y[46] + (y[47] + y[77] * (1.0 - parms[252])) * y[58] + (y[49] + y[75] * (1.0 - parms[252])) * y[58] + (y[50] + y[76] * (1.0 - parms[252])) * y[58] + y[51])/x[31])));
x[169] = -(y[55] * y[58] + y[56] * y[58] - y[23] * y[58] - y[16] * y[58] - (y[47] + y[77] * (1.0 - parms[252])) * y[58] - (y[49] + y[75] * (1.0 - parms[252])) * y[58] - (y[50] + y[76] * (1.0 - parms[252])) * y[58] - y[48] - y[52])/x[31];
x[170] = ((y[16] + y[23] + y[47] + y[77] * (1.0 - parms[252]) + y[49] + y[75] * (1.0 - parms[252]) + y[50] + y[76] * (1.0 - parms[252]) - y[56]) * y[58])/x[31];
x[173] = x[31]/(y[67] * y[58] * y[63]);
x[175] = (y[55] * y[58])/x[31];
x[176] = y[58] * (y[47] + y[49] + y[15] + y[16] + (y[77] + y[75] + y[76]) * (1.0 - parms[252]))/x[31];
x[177] = (y[17] + y[15] * y[58] + y[16] * y[58] + y[35])/x[31];
x[178] = (y[46] + (y[47] + y[49] + y[50] + (y[77] + y[75] + y[76]) * (1.0 - parms[252])) * y[58])/x[31];
x[207] = parms[59] * (parms[58] * y[12] * x[46] * (1.0 + parms[50]) - y[13]);
x[212] = 1.0/(1.0 + exp(-parms[64] * (x[171] - parms[67]))) * (parms[65] * (1.0 + parms[261] * parms[260] * parms[258] * (t > 4.0) * exp(parms[256] * (4.0 - t))) - parms[66]) + parms[66];
x[222] = 1.0/(1.0 + exp(-parms[106] * (x[171] - parms[109]))) * (parms[107] * (1.0 + parms[261] * parms[260] * parms[259] * (t > 4.0) * exp(parms[256] * (4.0 - t))) - parms[108]) + parms[108];
x[228] = parms[95] * (x[89] - y[28]);
ydot[13] = x[207];
ydot[28] = x[228];
ydot[73] = 2.0 * (x[15] - y[73]);
x[22] = std::max(0.064, 0.064 * (1.0 + parms[233] * tanh(parms[231] * (parms[232] - x[175]))));
x[48] = 1.0 - x[45]/y[68];
x[50] = parms[43] * x[49];
x[51] = parms[44] * x[49];
x[52] = parms[45] * (x[49]);
x[73] = (1.0 - parms[167]) * x[142];
x[93] = y[70] + x[92];
x[95] = y[70] + parms[103] * x[92];
x[104] = x[103] + parms[118] * x[102];
x[134] = parms[163] * x[31] * (1.0 - parms[239] * tanh(parms[237] * (parms[238] - x[175])));
x[143] = parms[167] * x[142] * (1.0 - parms[242] * tanh(parms[241] * (parms[240] - x[175])));
x[152] = y[70] + parms[182] * x[92];
x[201] = parms[41] * (x[35] - y[10]);
x[206] = (parms[46] * (x[203]/y[59]) + parms[47] * (x[45]/y[68] - parms[48]) + parms[49] * x[202]/y[11]) * y[12];
x[230] = std::max(parms[117] * x[17] * y[63] + y[57] * parms[230] * (1.0 + tanh(parms[229] * (parms[228] - x[175]))) - y[57], 0.0);
x[239] = parms[131] * x[102];
x[244] = ((parms[147] * (x[203]/y[59]) + parms[148] * x[202]/y[11]) * (1.0 - parms[236] * tanh(parms[234] * (parms[228] - x[175])))) * y[39];
x[254] = parms[181] * (x[151] - y[54]);
ydot[10] = x[201];
ydot[12] = x[206];
ydot[36] = x[239];
ydot[39] = x[244];
ydot[54] = x[254];
ydot[57] = x[230];
ydot[74] = parms[138] * (x[22] - y[74]);
x[74] = (1.0 - parms[146]) * x[104];
x[94] = x[93] * (1.0 + parms[104] * y[28]);
x[110] = parms[126] - parms[127] * x[90] - parms[128] * x[48];
x[114] = parms[132] + parms[133] * x[48];
x[125] = parms[146] * x[104];
x[127] = x[132] + (1.0 + parms[149]) * y[39] * x[135] + x[37] * x[136] + y[40] + x[13] * x[40] + x[143];
x[153] = (1.0 - parms[183]) * x[152];
x[241] = parms[164] * (x[134] - y[41]);
ydot[41] = x[241];
x[53] = x[49] - x[124] - x[87] * y[17] - x[94] * y[15] * y[58] - x[95] * y[16] * y[58] + x[85] * (1.0 - parms[90]) * y[13] - (x[52] + x[50] + x[51]);
x[76] = x[90] * y[35] + x[87] * y[17] + x[94] * y[15] * y[58] + x[150] * y[25] + y[53] * x[155] - (x[85] * y[43] + x[85] * (1.0 - parms[90]) * y[13] + x[85] * (1.0 - parms[91]) * y[34]) - x[93] * y[23] * y[58] + x[98] * y[21] * y[58] - y[30] * y[24] - x[37] * x[70] - parms[144] * x[63] - (1.0 + parms[81]) * y[20] * y[18] + x[74] - x[73] + x[67] + x[64];
x[113] = x[114] * y[67];
x[154] = (1.0 - parms[211]) * x[153];
x[54] = (1.0 - parms[135]) * x[53];
x[77] = (1.0 - parms[136]) * x[76];
x[115] = parms[134] * (x[113] - y[68]);
x[118] = parms[137] * x[102] + parms[135] * x[53] + parms[136] * x[76];
x[128] = x[150] * y[46] + x[152] * y[47] * y[58] + x[153] * y[49] * y[58] + y[50] * x[154] * y[58] + (y[53] + y[52]) * x[155] + (y[77] * x[152] + y[76] * x[154] + y[75] * x[153]) * (1.0 - parms[251]) * (1.0 - parms[252]) * y[58];
x[129] = x[150] * y[46] + x[152] * y[47] * y[58] + x[153] * y[49] * y[58] + y[50] * x[154] * y[58] + (y[53] + y[52]) * x[155] + (y[77] * x[152] + y[76] * x[154] + y[75] * x[153]) * y[58];
x[130] = x[152] * y[47] * y[58] + x[153] * y[49] * y[58] + y[50] * x[154] * y[58] + (y[77] * x[152] + y[76] * x[154] + y[75] * x[153]) * (1.0 - parms[251]) * (1.0 - parms[252]) * y[58];
x[131] = x[152] * y[47] * y[58] + x[153] * y[49] * y[58] + y[50] * x[154] * y[58] + (y[77] * x[152] + y[76] * x[154] + y[75] * x[153]) * y[58];
x[187] = std::min(1.0, x[13] * x[40]/(y[58] * (x[152] * y[47] + x[153] * y[49] + y[50] * x[154]) * (1.0 - parms[251]) * (1.0 - parms[252]))) * parms[260];
ydot[68] = x[115];
x[55] = x[54] - x[208] * y[58] - x[207];
x[56] = x[54]/(x[38] * y[2]);
x[80] = x[77] - x[79] - x[183];
x[117] = x[118] + x[119] + x[120] + x[121] + x[122] + x[123];
x[126] = x[127] + x[128];
x[181] = (x[87] * y[17] + x[94] * y[15] * y[58] + x[95] * y[16] * y[58])/x[54];
x[188] = x[187] * y[49];
x[189] = x[187] * y[50];
x[190] = x[187] * y[47];
ydot[75] = parms[253] * (x[188] - y[75]);
ydot[76] = parms[255] * (x[189] - y[76]);
ydot[77] = parms[254] * (x[190] - y[77]);
x[11] = (x[12] + parms[8] * (x[56] - x[202]/y[11])) * y[2];
x[57] = (1.0 - parms[57]) * x[55];
x[61] = parms[57] * x[55] - x[182];
x[81] = x[80];
x[147] = parms[168] * x[126];
x[148] = parms[170] * x[126];
x[255] = parms[253] * (x[188] - y[75]);
x[256] = parms[255] * (x[189] - y[76]);
x[257] = parms[254] * (x[190] - y[77]);
x[58] = std::max(0.0, (parms[53] + parms[54] * (y[72] * y[65])/(x[31]/y[58])) * x[57]);
x[59] = std::max(0.0, (parms[55] + parms[56] * (y[72] * y[65])/(x[31]/y[58])) * x[57]);
x[62] = x[38] * y[4] - x[61];
x[193] = parms[9] * (x[11] - y[4]);
x[245] = parms[169] * (x[147] - y[43]);
x[246] = parms[171] * (x[148] - y[44]);
x[249] = -parms[177] * (y[71] * x[157]/y[58]) - x[257];
x[250] = -(1.0 - parms[177]) * (y[71] * x[157]/y[58]) - x[255];
x[251] = parms[208] * parms[209] * x[13] * x[40]/y[58] - x[256];
ydot[4] = x[193];
ydot[43] = x[245];
ydot[44] = x[246];
ydot[47] = x[249];
ydot[49] = x[250];
ydot[50] = x[251];
x[60] = x[57] - x[58] - x[59];
x[116] = x[117] + x[124] + x[51] + x[125] + x[85] * y[43] + x[100] * y[44] + x[59] - x[184] + x[99];
x[158] = (x[160] + x[97] * y[57] + x[98] * y[21] + x[159] - x[154] * y[50] - x[152] * y[47] - x[153] * y[49] - x[93] * y[23] - x[95] * y[16] - (y[77] * x[152] + y[76] * x[154] + y[75] * x[153]) * (1.0 - parms[251]) * (1.0 - parms[252])) * y[58] - x[58] - x[82] + x[186] - x[150] * y[48] - y[52] * x[155];
x[209] = parms[62] * (x[62]/y[58]);
x[211] = (1.0 - x[212]) * parms[63] * (x[62]/y[58]);
ydot[16] = x[211];
x[101] = (1.0 - parms[137]) * x[102] + x[52] + x[103] + x[143] + x[73] - x[104] - x[90] * y[35] + x[85] * (1.0 - parms[91]) * y[34] + x[60] + x[81] + x[160] * y[58] + x[50] - x[68] - x[65] + x[185];
x[144] = x[126] - x[116] - x[133];
x[156] = -(x[157] + x[158])/x[31];
x[210] = (1.0 - x[222]) * x[209];
x[220] = (parms[105] * x[226]/y[58] + x[209]);
ydot[15] = x[210];
x[91] = parms[96] + parms[97]/(1.0 + exp(-parms[98] * (y[35]/x[101])));
x[108] = y[38] * x[101];
x[109] = x[110] * x[101];
x[111] = x[101] - y[31];
x[138] = std::max(parms[156], parms[155] - parms[157] * x[139] * (x[144] - x[145])/(y[42] * x[38]) * parms[262]);
x[146] = x[144] + x[245] + x[246] + x[247] * y[58];
x[149] = parms[174] + parms[175] * x[156];
x[179] = x[144]/x[31];
x[213] = x[62] - x[210] * y[58] - x[211] * y[58] - x[165];
x[217] = x[220] - x[219] - x[210];
x[221] = (1.0 - x[222]) * x[220];
ydot[17] = x[213];
ydot[23] = x[221];
x[112] = y[33] - x[111];
x[140] = x[138] * y[42];
x[167] = x[17] * y[63] + x[152] * y[47] + x[153] * y[49] + x[93] * y[23] + x[95] * y[16] + x[154] * y[50] + x[58]/y[58] + x[82]/y[58] + x[217] + x[261] + x[150] * y[48]/y[58] + y[52] * x[155]/y[58] + (y[77] * x[152] + y[76] * x[154] + y[75] * x[153]) * (1.0 - parms[251]) * (1.0 - parms[252]);
x[168] = x[26]/y[58] + x[160] + x[186]/y[58] + x[97] * y[57] + x[98] * y[21] + x[161]/y[58] + x[249] + x[257] + x[250] + x[255] + x[211] + x[221] + x[251] + x[256] + x[159] + x[258]/y[58] + x[259]/y[58] - x[230];
x[180] = (x[90] * y[35])/x[111];
x[229] = parms[99] * (x[91] - y[29]);
x[233] = parms[129] * (x[109] - y[33]);
x[237] = parms[124] * (x[108] - y[69]);
x[248] = parms[176] * (x[149] - y[71]);
x[253] = x[146] - x[249] * y[58] - x[250] * y[58] - x[251] * y[58] - x[255] * y[58] - x[256] * y[58] - x[257] * y[58] - x[252] - parms[205] * x[13] * x[40]/y[58];
x[260] = (x[157]/y[58] + x[158]/y[58] + x[161]/y[58] + x[249] + x[258]/y[58] + x[259]/y[58] + x[250] + x[211] + x[221] + x[251] + x[255] + x[257] + x[256] + x[159] - x[261]);
ydot[29] = x[229];
ydot[33] = x[233];
ydot[46] = x[253];
ydot[55] = x[260];
ydot[69] = x[237];
ydot[71] = x[248];
x[105] = x[106] * x[101] + x[107] * (y[34] + y[36]) + x[237];
x[141] = x[140] * x[38] + parms[250] * x[31];
x[218] = x[260] - x[230];
x[223] = x[253] - x[258];
x[235] = x[237] + x[236];
x[262] = parms[194] * ((x[167] - x[168])/x[168]);
ydot[21] = x[218];
ydot[25] = x[223];
ydot[35] = x[235];
ydot[58] = x[262];
x[232] = parms[119] * (x[105] - y[31]);
x[240] = x[111] - y[33] + x[235] - x[239];
x[242] = parms[162] * (x[141] - y[40]);
ydot[31] = x[232];
ydot[34] = x[240];
ydot[40] = x[242];
x[84] = (x[213] + x[235] + x[223] + x[224]) + x[72] + parms[110] * (x[245] + x[240] + x[207]) - (x[245] + x[240] + x[207] + x[226] + x[166]) - x[239] + (x[219] * y[58] + x[210] * y[58] - x[221] * y[58] + x[218] * y[58]);
x[225] = parms[110] * ((x[245] + x[240] + x[207]));
ydot[26] = x[225];
x[227] = std::max(x[84], -y[24]);
ydot[24] = x[227];
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
parms[220] = 0.025; 
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
