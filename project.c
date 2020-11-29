// Gavin Bowman
// Codon Bias Measures Research Project
//
// This program reads a gene sequence and reference
// table from txt and csv files respectively. It can then
// display RSCU values for each codon, as well as the 
// CAI and SCUO values for the whole sequence. Future
// goals include the ability to parse multiple sequences
// and more Codon Bias Measures.

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>

// temporary struct for parsing data from csv
// data will be moved to array of amino structs
struct myTemps {
	char name[3];
	char trip[3];
	double freq;
};

// struct for each codon of a amino acid
struct Amino {
	char* name; // name of amino acid
	char* trip; // nucleotide sequence for this codon
	double freq; //frequecny per thousand from reference set
	int count; // number of occurences in this organism
	double rscu; //calculated rscu value
	double pij; // intermediate value for scuo calculation
	double hij; // ditto
	double hi; // ""
	double hiMax; // ""
};

// counts occurences of each codon of given amino acid across sequence
// used in for loop; bool return value tells loop to continue
bool countCodon(struct Amino** list, char* seq, int size) {

	int i = 0;
	for(; i < size; i++) {
		if(strcmp(seq, list[i]->trip) == 0) {
			list[i]->count++;
			return true;
			break;
		}
		else {
		}
	}
	return false;
}

// calculates the rscu values for all codons of a given amino acid
// also calculates pij and hij values which will later be important
// for scuo.
void rscu(struct Amino** list, int size) {

	int i = 0;
	double tempList[size]; // for easier access to count of each codon

	for(; i < size; i++) {
		tempList[i] = (double)(list[i]->count);
	}

	// for each codon of this amino acid...
	for(i = 0; i < size; i++) {

		double rscu = 0;
		double tempSum = tempList[i]; // initialize sum as this codon's count
		int j = 0;

		// for each codon of this amino acid, we add the counts
		// taking into account that we already know that of the 
		// current one (i).
		for(; j < size; j++) {
			if(i == j) {
				continue;
			}
			else {
				tempSum = tempSum + tempList[j];
			}
		}

		// since this sum is important for scuo, we will
		// use what we need of it here.
		//
		// we calculate the pij and hij values for each
		// codon of the amino acid using the sum we just 
		// got for our rscu calculation.
		if(tempList[i] && tempSum != 0) {
			list[i]->pij = tempList[i] / tempSum;
               		list[i]->hij = (-1 * list[i]->pij) * log10(list[i]->pij);
		}
		else {
			// note - I put this separate because for
			// whatever reason, log10(0) is not pretty
			list[i]->pij = 0;
			list[i]->hij = 0;
		}
		
		double denom;
		if(tempSum != 0) {
			denom = (1/(double)(size)) * tempSum;
		}
		else {
			continue;
		}
		rscu = (tempList[i]) / denom;
		list[i]->rscu = rscu;
	}				
}

// misnamed. should be fiNumerator. nevertheless,
// this function calculates the numerator for the 
// fi value for each amino acid. this is an intermediate
// value for calculating scuo.
double fiDenom(struct Amino** list, int size) {

	int i = 1;
	double sum = list[0]->count;

	for(; i < size; i++) {
		sum = sum + list[i]->count;
	}
	return sum;
}

// calculates the hi value for each amino acid. this
// is called below in the calcOi function.
double calcHi(struct Amino** list, int size) {

	int i = 0;
	double sum = 0;

	// sum of absolute value of hij values for all
	// codons of this amino acid. originally tried to
	// use the fabs() fucntion on the hij value for each
	// codon but was getting weird values (inf, nan, and 
	// numbers that were just objectively incorrect). as
	// such, I ended up going the more verbose route which
	// resulted in better calculations.
	for(; i < size; i++) {
		if(list[i]->pij != 0) {
			sum = sum + (list[i]->pij * log10(list[i]->pij));
		}
		else sum = sum + 0; // similar to rscu calculation, had to
				    // do this rather than take log of 0.
	}
	return (sum * -1);
}

// calculates the oi value for each amino acid. I had originally
// tried to calculate hiMax as (-1 * log10(1 / size)), but the
// compiler didn't like there to be division inside the log function.
// as such, I hardcoded the hiMax for each amino acid based on its
// number of synonymous codons, and put them in arrays with their
// index correlating to said number.
double calcOi(double hi, int size) {

	double hiMax[7] = {0, 0, (-1 * log10(0.5)), (-1 * log10(0.333333)), (-1 * log10(0.25)), 0, (-1 * log10(0.16666666)) };
	double retval = ((hiMax[size] - hi) / hiMax[size]);
	return retval;
}

double calcFiOi(struct Amino** list, int size) {

	double hi = calcHi(list, size);
	double oi = calcOi(hi, size);
}

// displays the rscu value for each codon
void dispRscu(struct Amino** gly, struct Amino** glu, struct Amino** asp, struct Amino** val, struct Amino** ala, struct Amino** arg, struct Amino** ser, struct Amino** lys, struct Amino** asn, struct Amino** met, struct Amino** ile, struct Amino** thr, struct Amino** trp, struct Amino** en, struct Amino** cys, struct Amino** tyr, struct Amino** leu, struct Amino** phe, struct Amino** gln, struct Amino** his, struct Amino** pro) {

	// for each amino acid, print each codon's nucleotide sequence and 
	// its rscu value.

	int i = 0;
	for(; i < 4; i++) {
		printf("%s's RSCU value is %f    ", gly[i]->trip, gly[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 2; i++) {
		printf("%s's RSCU value is %f    ", glu[i]->trip, glu[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 2; i++) {
		printf("%s's RSCU value is %f    ", asp[i]->trip, asp[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 4; i++) {
		printf("%s's RSCU value is %f    ", val[i]->trip, val[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 4; i++) {
		printf("%s's RSCU value is %f    ", ala[i]->trip, ala[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 6; i++) {
		if(i == 4) {
			printf("\n");
		}
		printf("%s's RSCU value is %f    ", arg[i]->trip, arg[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 6; i++) {
		if(i == 4) {
			printf("\n");
		}
		printf("%s's RSCU value is %f    ", ser[i]->trip, ser[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 2; i++) {
		printf("%s's RSCU value is %f    ", lys[i]->trip, lys[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 2; i++) {
		printf("%s's RSCU value is %f    ", asn[i]->trip, asn[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 1; i++) {
		printf("%s's RSCU value is %f    ", met[i]->trip, met[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 3; i++) {
		printf("%s's RSCU value is %f    ", ile[i]->trip, ile[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 4; i++) {
		printf("%s's RSCU value is %f    ", thr[i]->trip, thr[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 1; i++) {
		printf("%s's RSCU value is %f    ", trp[i]->trip, trp[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 3; i++) {
		printf("%s's RSCU value is %f    ", en[i]->trip, en[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 2; i++) {
		printf("%s's RSCU value is %f    ", cys[i]->trip, cys[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 2; i++) {
		printf("%s's RSCU value is %f    ", tyr[i]->trip, tyr[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 6; i++) {
		if(i == 4) {
			printf("\n");
		}
		printf("%s's RSCU value is %f    ", leu[i]->trip, leu[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 2; i++) {
		printf("%s's RSCU value is %f    ", phe[i]->trip, phe[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 2; i++) {
		printf("%s's RSCU value is %f    ", gln[i]->trip, gln[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 2; i++) {
		printf("%s's RSCU value is %f    ", his[i]->trip, his[i]->rscu);
	}
	printf("\n\n");
	for(i = 0; i < 4; i++) {
		printf("%s's RSCU value is %f    ", pro[i]->trip, pro[i]->rscu);
	}
	printf("\n\n");
}

// returns a number from 0 to 20, indicating which
// amino acid the given codon codes for.
// this is used when parsing the gene sequence.
int incIndex(char* other) {
		if(strcmp(other, "GGG") == 0 || strcmp(other, "GGA") == 0 || strcmp(other, "GGT") == 0 || strcmp(other, "GGC") == 0) {
			return 0;
		}
		else if(strcmp(other, "GAG") == 0 || strcmp(other, "GAA") == 0) {
			return 1;
		}
		else if(strcmp(other, "GAT") == 0 || strcmp(other, "GAC") == 0|| strcmp(other, "GAU") == 0) {
			return 2;
		}
		else if(strcmp(other, "GTG") == 0 || strcmp(other, "GTA") == 0 || strcmp(other, "GTT") == 0 || strcmp(other, "GTC") == 0 || strcmp(other, "GUU") == 0 || strcmp(other, "GUC") == 0 || strcmp(other, "GUA") == 0 || strcmp(other, "GUG") == 0) {
			return 3;
		}
		else if(strcmp(other, "GCG") == 0 || strcmp(other, "GCA") == 0 || strcmp(other, "GCT") == 0 || strcmp(other, "GCC") == 0) {
			return 4;
		}
		else if(strcmp(other, "AGG") == 0 || strcmp(other, "AGA") == 0 || strcmp(other, "CGG") == 0 || strcmp(other, "CGA") == 0 || strcmp(other, "CGT") == 0 || strcmp(other, "CGC") == 0 || strcmp(other, "CGU") == 0 || strcmp(other, "CGC") == 0 || strcmp(other, "CGA") == 0 || strcmp(other, "CGG") == 0) {
			return 5;
		}
		else if(strcmp(other, "AGT") == 0 || strcmp(other, "AGC") == 0 || strcmp(other, "TCG") == 0 || strcmp(other, "TCA") == 0 || strcmp(other, "TCT") == 0 || strcmp(other, "TCC") == 0 || strcmp(other, "UCU") == 0 || strcmp(other, "UCC") == 0 || strcmp(other, "UCA") == 0 || strcmp(other, "UCG") == 0 || strcmp(other, "AGU") == 0 || strcmp(other, "AGC") == 0) {
			return 6;
		}
		else if(strcmp(other, "AAG") == 0 || strcmp(other, "AAA") == 0) {
			return 7;
		}
		else if(strcmp(other, "AAT") == 0 || strcmp(other, "AAC") == 0 || strcmp(other, "AAU") == 0) {
			return 8;
		}
		else if(strcmp(other, "ATG") == 0 || strcmp(other, "AUG") == 0) {
			return 9;
		}
		else if(strcmp(other, "ATA") == 0 || strcmp(other, "ATT") == 0 || strcmp(other, "ATC") == 0 || strcmp(other, "AUU") == 0 || strcmp(other, "AUC") == 0 || strcmp(other, "AUA") == 0) {
			return 10;
		}
		else if(strcmp(other, "ACG") == 0 || strcmp(other, "ACA") == 0 || strcmp(other, "ACT") == 0 || strcmp(other, "ACC") == 0 || strcmp(other, "ACU") == 0) {
			return 11;
		}
		else if(strcmp(other, "TGG") == 0 || strcmp(other, "UGG") == 0) {
			return 12;
		}
		else if(strcmp(other, "TGA") == 0 || strcmp(other, "TAG") == 0 || strcmp(other, "TAA") == 0 || strcmp(other, "UAA") == 0 || strcmp(other, "UGA") == 0 || strcmp(other, "UAG") == 0) {
			return 13;
		}
		else if(strcmp(other, "TGT") == 0 || strcmp(other, "TGC") == 0 || strcmp(other, "UGU") == 0 || strcmp(other, "UGC") == 0) {
			return 14;
		}
		else if(strcmp(other, "TAT") == 0 || strcmp(other, "TAC") == 0 || strcmp(other, "UAU") == 0 || strcmp(other, "UAC") == 0) {
			return 15;
		}
		else if(strcmp(other, "TTG") == 0 || strcmp(other, "TTA") == 0 || strcmp(other, "CTG") == 0 || strcmp(other, "CTA") == 0 || strcmp(other, "CTT") == 0 || strcmp(other, "CTC") == 0 || strcmp(other, "UUA") == 0 || strcmp(other, "UUG") == 0 || strcmp(other, "CUU") == 0 || strcmp(other, "CUC") == 0 || strcmp(other, "CUA") == 0 || strcmp(other, "CUG") == 0) {
			return 16;
		}
		else if(strcmp(other, "TTT") == 0 || strcmp(other, "TTC") == 0 || strcmp(other, "UUU") == 0 || strcmp(other, "UUC") == 0) {
			return 17;
		}
		else if(strcmp(other, "CAG") == 0 || strcmp(other, "CAA") == 0) {
			return 18;
		}
		else if(strcmp(other, "CAT") == 0 || strcmp(other, "CAC") == 0 || strcmp(other, "CAU") == 0) {
			return 19;
		}
		else if(strcmp(other, "CCG") == 0 || strcmp(other, "CCA") == 0 || strcmp(other, "CCT") == 0 || strcmp(other, "CCC") == 0 || strcmp(other, "CCU") == 0) {
			return 20;
		}
		else {
		}

}

// calculates the cai value for the organism in question.
double cai(struct Amino** gly, struct Amino** glu, struct Amino** asp, struct Amino** val, struct Amino** ala, struct Amino** arg, struct Amino** ser, struct Amino** lys, struct Amino** asn, struct Amino** met, struct Amino** ile, struct Amino** thr, struct Amino** trp, struct Amino** en, struct Amino** cys, struct Amino** tyr, struct Amino** leu, struct Amino** phe, struct Amino** gln, struct Amino** his, struct Amino** pro, int len, char* seq) {

        double glyM, gluM, aspM, valM, alaM, argM, serM, lysM, asnM, metM, ileM, thrM, trpM, enM, cysM, tyrM, leuM, pheM, glnM, hisM, proM = 0;

	double cai = 0;
        double values[len];
        int i = 0;

	// for each amino acid in the sequnce, the frequency of the most frequent
	// codon is determined, and stored in the values array.
        for(; i < 4; i++) {
                if(gly[i]->freq > glyM) {
                        glyM = gly[i]->freq;
                }
        }
        for(i = 0; i < 2; i++) {
                if(glu[i]->freq > gluM) {
                        gluM = glu[i]->freq;
                }

        }
        for(i = 0; i < 2; i++) {
                if(asp[i]->freq > aspM) {
                        aspM = asp[i]->freq;
                }
        }
        for(i = 0; i < 4; i++) {
                if(val[i]->freq > valM) {
                        valM = val[i]->freq;
                }
        }
        for(i = 0; i < 4; i++) {
                if(ala[i]->freq > alaM) {
                        alaM = ala[i]->freq;
                }
        }
        for(i = 0; i < 6; i++) {
                if(arg[i]->freq > argM) {
                        argM = arg[i]->freq;
                }
        }
        for(i = 0; i < 6; i++) {
                if(ser[i]->freq > serM) {
                        serM = ser[i]->freq;
                }
        }
        for(i = 0; i < 2; i++) {
                if(lys[i]->freq > lysM) {
                        lysM = lys[i]->freq;
                }
        }
        for(i = 0; i < 2; i++) {
                if(asn[i]->freq > asnM) {
                        asnM = asn[i]->freq;
                }
        }
        for(i = 0; i < 1; i++) {
                if(met[i]->freq > metM) {
                        metM = met[i]->freq;
                }
        }
        for(i = 0; i < 3; i++) {
                if(ile[i]->freq > ileM) {
                        ileM = ile[i]->freq;
                }
        }
        for(i = 0; i < 4; i++) {
                if(thr[i]->freq > thrM) {
                        thrM = thr[i]->freq;
                }
        }
        for(i = 0; i < 1; i++) {
                if(trp[i]->freq > trpM) {
                        trpM = trp[i]->freq;
                }
        }
        for(i = 0; i < 3; i++) {
                if(en[i]->freq > enM) {
                        enM = en[i]->freq;
                }
        }
        for(i = 0; i < 2; i++) {
                if(cys[i]->freq > cysM) {
                        cysM = cys[i]->freq;
                }
        }
        for(i = 0; i < 2; i++) {
                if(tyr[i]->freq > tyrM) {
                        tyrM = tyr[i]->freq;
                }
        }
        for(i = 0; i < 6; i++) {
                if(leu[i]->freq > leuM) {
                        leuM = leu[i]->freq;
                }
        }
        for(i = 0; i < 2; i++) {
                if(phe[i]->freq > pheM) {
                        pheM = phe[i]->freq;
                }
        }
        for(i = 0; i < 2; i++) {
                if(gln[i]->freq > glnM) {
                        glnM = gln[i]->freq;
                }
	}
        for(i = 0; i < 2; i++) {
                if(his[i]->freq > hisM) {
                        hisM = his[i]->freq;
                }
        }
        for(i = 0; i < 4; i++) {
                if(pro[i]->freq > proM) {
                        proM = pro[i]->freq;
                }
        }
        
        char sequence[strlen(seq)];
        strcpy(sequence, seq);
        char other[4];
        char rep[4] = "";
        int loop = 1;
        int g = 0;

	// the frequency of each codon in the sequence is divided by
	// that of its most frequent synonymous codon.
        while(((strlen(sequence) - 1) % 3) == 0 && strlen(sequence) != 1) {
                strncpy(other, sequence, 4);
                other[3] = '\0';
                if(strcmp(other, "") == 0) {
                        break;
                }
                else {
                        if(incIndex(other) == 0) {
                                for(g = 0; g < 4; g++) {
                                        if(strcmp(other, gly[g]->trip) == 0) {
                                                values[loop - 1] = (gly[g]->freq) / glyM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 1) {
				for(g = 0; g < 2; g++) {
					if(strcmp(other, glu[g]->trip) == 0) {
						values[loop - 1] = (glu[g]->freq) / gluM;
						strcpy(sequence, rep);
						strncpy(sequence, &seq[(3 * loop)], strlen(seq));
						loop++;
						strcpy(other, rep);
						break;
					}
				}
			}
			else if(incIndex(other) == 2) {
                                for(g = 0; g < 2; g++) {
                                        if(strcmp(other, asp[g]->trip) == 0) {
                                                values[loop - 1] = (asp[g]->freq) / aspM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 3) {
                                for(g = 0; g < 4; g++) {
                                        if(strcmp(other, val[g]->trip) == 0) {
                                                values[loop - 1] = (val[g]->freq) / valM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 4) {
                                for(g = 0; g < 4; g++) {
                                        if(strcmp(other, ala[g]->trip) == 0) {
                                                values[loop - 1] = (ala[g]->freq) / alaM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 5) {
                                for(g = 0; g < 6; g++) {
                                        if(strcmp(other, arg[g]->trip) == 0) {
                                                values[loop - 1] = (arg[g]->freq) / alaM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 6) {
                                for(g = 0; g < 6; g++) {
                                        if(strcmp(other, ser[g]->trip) == 0) {
                                                values[loop - 1] = (ser[g]->freq) / serM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 7) {
                                for(g = 0; g < 2; g++) {
                                        if(strcmp(other, lys[g]->trip) == 0) {
                                                values[loop - 1] = (lys[g]->freq) / lysM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 8) {
                                for(g = 0; g < 2; g++) {
                                        if(strcmp(other, asn[g]->trip) == 0) {
                                                values[loop - 1] = (asn[g]->freq) / asnM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 9) {
                                for(g = 0; g < 1; g++) {
                                        if(strcmp(other, met[g]->trip) == 0) {
                                                values[loop - 1] = (met[g]->freq) / metM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 10) {
                                for(g = 0; g < 3; g++) {
                                        if(strcmp(other, ile[g]->trip) == 0) {
                                                values[loop - 1] = (ile[g]->freq) / ileM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 11) {
                                for(g = 0; g < 4; g++) {
                                        if(strcmp(other, thr[g]->trip) == 0) {
                                                values[loop - 1] = (thr[g]->freq) / thrM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 12) {
                                for(g = 0; g < 1; g++) {
                                        if(strcmp(other, trp[g]->trip) == 0) {
                                                values[loop - 1] = (trp[g]->freq) / trpM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 13) {
                                for(g = 0; g < 3; g++) {
                                        if(strcmp(other, en[g]->trip) == 0) {
                                                values[loop - 1] = (en[g]->freq) / enM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 14) {
                                for(g = 0; g < 2; g++) {
                                        if(strcmp(other, cys[g]->trip) == 0) {
                                                values[loop - 1] = (cys[g]->freq) / cysM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 15) {
                                for(g = 0; g < 2; g++) {
                                        if(strcmp(other, tyr[g]->trip) == 0) {
                                                values[loop - 1] = (tyr[g]->freq) / tyrM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 16) {
                                for(g = 0; g < 6; g++) {
                                        if(strcmp(other, leu[g]->trip) == 0) {
                                                values[loop - 1] = (leu[g]->freq) / leuM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 17) {
                                for(g = 0; g < 2; g++) {
                                        if(strcmp(other, phe[g]->trip) == 0) {
                                                values[loop - 1] = (phe[g]->freq) / pheM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 18) {
                                for(g = 0; g < 2; g++) {
                                        if(strcmp(other, gln[g]->trip) == 0) {
                                                values[loop - 1] = (gln[g]->freq) / glnM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 19) {
                                for(g = 0; g < 2; g++) {
                                        if(strcmp(other, his[g]->trip) == 0) {
                                                values[loop - 1] = (his[g]->freq) / hisM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
			else if(incIndex(other) == 20) {
                                for(g = 0; g < 4; g++) {
                                        if(strcmp(other, pro[g]->trip) == 0) {
                                                values[loop - 1] = (pro[g]->freq) / proM;
                                                strcpy(sequence, rep);
                                                strncpy(sequence, &seq[(3 * loop)], strlen(seq));
                                                loop++;
                                                strcpy(other, rep);
                                                break;
                                        }
                                }
                        }
                }
	}

	// final cai value is determined by taking the geometric
	// mean of the values calculated by the previous loop.
	int qwer = 1;
	cai = sqrt(values[0]);
	for(; qwer < len; qwer++) {
		if(values[qwer] != 0) {

			// multiplying the values as they stood
			// produced results too small for a double
			// as such, I opted to square them here,
			// then halve the denominator in the power.
			cai = cai * sqrt(values[qwer]);
		}
	}

	// in other words, normally cai is calculated as
	// (val1*val2*...*valn)^(1/n)
	// to prevent arithmetic underflow, I instead implemented
	// it as 
	// (val1^2*val2^2*...*valn^2)^(1/(n/2))
	// which produces the same final result
	cai = pow(cai, (1 / ((double)(len) / 2)));
	return cai;
}

int main() {

	char buffer[800];

	// get names of input files
	printf("\nWhat is the name of the file containing the sequence (.txt)?\n\n > ");

	char in1[20];
	scanf("%s", in1);
	printf("\n");

	FILE *file = fopen(in1, "r");
	
	if(file == NULL) {
		fputs("Failed to open the file. Please check that you entered the right file name.\n", stderr);
	}
	while (fgets(buffer, sizeof(buffer), file) != NULL) {
	}
	fclose(file);

	printf("What is the name of the file containing the reference organism for CAI calculation (.csv)?\n\n > ");

	char in2[20];
	scanf("%s", in2);
	printf("\n");

	FILE* fh = fopen(in2, "r");

	if(fh == NULL) {
		fputs("Failed to open the file. Please check that you entered the right file name.\n", stderr);
	}

	struct myTemps temps[100];
	size_t count = 0;

	// read raw data from csv into unorganized structs containing
	// the name, nucleotide sequence and frequency separated into
	// their own variables.
	for(; count < sizeof(temps)/sizeof(temps[0]); ++count) {
		int got = fscanf(fh, "%[^,],%[^,],%lf ", temps[count].name, temps[count].trip, &temps[count].freq);
		if(got != 3) break;
	}
	fclose(fh);

	int counter[21] = { 0 };
	char sequence[strlen(buffer)];
	strcpy(sequence, buffer);
	char other[4];
	char rep[4] = "";
	int loop = 1;

	// count the occurences of each amino acid across the sequence
	while(((strlen(sequence) - 1) % 3) == 0) {
		strncpy(other, sequence, 4);
		other[3] = '\0';
		if(strcmp(other, "") == 0) {
			break;
		}
		else{
			counter[incIndex(other)]++;
		}
		strcpy(sequence, rep);
		strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
		loop++;
		strcpy(other, rep);
	}

	// initialize arrays of codons for each amino acid
	struct Amino* arrGly[4];
	struct Amino* arrGlu[2];
	struct Amino* arrAsp[2];
	struct Amino* arrVal[4];
	struct Amino* arrAla[4];
	struct Amino* arrArg[6];
	struct Amino* arrSer[6];
	struct Amino* arrLys[2];
	struct Amino* arrAsn[2];
	struct Amino* arrMet[1];
	struct Amino* arrIle[3];
	struct Amino* arrThr[4];
	struct Amino* arrTrp[1];
	struct Amino* arrEnd[3];
	struct Amino* arrCys[2];
	struct Amino* arrTyr[2];
	struct Amino* arrLeu[6];
	struct Amino* arrPhe[2];
	struct Amino* arrGln[2];
	struct Amino* arrHis[2];
	struct Amino* arrPro[4];

	// counts for each amino acid.
	int glyC, gluC, aspC, valC, alaC, argC, serC, lysC, asnC, metC, ileC, thrC, trpC, endC, cysC, tyrC, leuC, pheC, glnC, hisC, proC = 0;

	int incr = 0;
	char tempBuff[4];

	// populate the arrays with name, nucleotide sequence, and frequency for
	// each codon.
	for(; incr < 64; incr++) {
		strncpy(tempBuff, temps[incr].name, 4);
		tempBuff[3] = '\0';
		
		if(strcmp(tempBuff, "Gly") == 0) {
			struct Amino *aminoGly = (struct Amino*) malloc(sizeof(struct Amino));
			aminoGly->name = "Gly";
			aminoGly->trip = temps[incr].trip;
			aminoGly->freq = temps[incr].freq;
			aminoGly->rscu = 0;
			aminoGly->count = 0;
			arrGly[glyC] = aminoGly;
			glyC++;
		}
		else if(strcmp(tempBuff, "Glu") == 0) {
			struct Amino *aminoGlu = (struct Amino*) malloc(sizeof(struct Amino));
			aminoGlu->name = "Glu";
			aminoGlu->trip = temps[incr].trip;
			aminoGlu->freq = temps[incr].freq;
			aminoGlu->rscu = 0;
			aminoGlu->count = 0;
			arrGlu[gluC] = aminoGlu;
			gluC++;
		}
		else if(strcmp(tempBuff, "Asp") == 0) {
			struct Amino *aminoAsp = (struct Amino*) malloc(sizeof(struct Amino));
			aminoAsp->name = "Asp";
			aminoAsp->trip = temps[incr].trip;
			aminoAsp->freq = temps[incr].freq;
			aminoAsp->rscu = 0;
			aminoAsp->count = 0;
			arrAsp[aspC] = aminoAsp;
			aspC++;
		}	
		else if(strcmp(tempBuff, "Val") == 0) {
			struct Amino *aminoVal = (struct Amino*) malloc(sizeof(struct Amino));
			aminoVal->name = "Val";
			aminoVal->trip = temps[incr].trip;
			aminoVal->freq = temps[incr].freq;
			aminoVal->rscu = 0;
			aminoVal->count = 0;
			arrVal[valC] = aminoVal;
			valC++;
		}
		else if(strcmp(tempBuff, "Ala") == 0) {
			struct Amino *aminoAla = (struct Amino*) malloc(sizeof(struct Amino));
			aminoAla->name = "Ala";
			aminoAla->trip = temps[incr].trip;
			aminoAla->freq = temps[incr].freq;
			aminoAla->rscu = 0;
			aminoAla->count = 0;
			arrAla[alaC] = aminoAla;
			alaC++;
		}
		else if(strcmp(tempBuff, "Arg") == 0) {
			struct Amino *aminoArg = (struct Amino*) malloc(sizeof(struct Amino));
			aminoArg->name = "Arg";
			aminoArg->trip = temps[incr].trip;
			aminoArg->freq = temps[incr].freq;
			aminoArg->rscu = 0;
			aminoArg->count = 0;
			arrArg[argC] = aminoArg;
			argC++;
		}
		else if(strcmp(tempBuff, "Ser") == 0) {
			struct Amino *aminoSer = (struct Amino*) malloc(sizeof(struct Amino));
			aminoSer->name = "Ser";
			aminoSer->trip = temps[incr].trip;
			aminoSer->freq = temps[incr].freq;
			aminoSer->rscu = 0;
			aminoSer->count = 0;
			arrSer[serC] = aminoSer;
			serC++;
		}
		else if(strcmp(tempBuff, "Lys") == 0) {
			struct Amino *aminoLys = (struct Amino*) malloc(sizeof(struct Amino));
			aminoLys->name = "Lys";
			aminoLys->trip = temps[incr].trip;
			aminoLys->freq = temps[incr].freq;
			aminoLys->rscu = 0;
			aminoLys->count = 0;
			arrLys[lysC] = aminoLys;
			lysC++;
		}
		else if(strcmp(tempBuff, "Asn") == 0) {
			struct Amino *aminoAsn = (struct Amino*) malloc(sizeof(struct Amino));
			aminoAsn->name = "Asn";
			aminoAsn->trip = temps[incr].trip;
			aminoAsn->freq = temps[incr].freq;
			aminoAsn->rscu = 0;
			aminoAsn->count = 0;
			arrAsn[asnC] = aminoAsn;
			asnC++;
		}
		else if(strcmp(tempBuff, "Met") == 0) {
			struct Amino *aminoMet = (struct Amino*) malloc(sizeof(struct Amino));
			aminoMet->name = "Met";
			aminoMet->trip = temps[incr].trip;
			aminoMet->freq = temps[incr].freq;
			aminoMet->rscu = 1;
			aminoMet->count = 0;
			arrMet[0] = aminoMet;
		}
		else if(strcmp(tempBuff, "Ile") == 0) {
			struct Amino *aminoIle = (struct Amino*) malloc(sizeof(struct Amino));
			aminoIle->name = "Ile";
			aminoIle->trip = temps[incr].trip;
			aminoIle->freq = temps[incr].freq;
			aminoIle->rscu = 0;
			aminoIle->count = 0;
			arrIle[ileC] = aminoIle;
			ileC++;
		}
		else if(strcmp(tempBuff, "Thr") == 0) {
			struct Amino *aminoThr = (struct Amino*) malloc(sizeof(struct Amino));
			aminoThr->name = "Thr";
			aminoThr->trip = temps[incr].trip;
			aminoThr->freq = temps[incr].freq;
			aminoThr->rscu = 0;
			aminoThr->count = 0;
			arrThr[thrC] = aminoThr;
			thrC++;
		}
		else if(strcmp(tempBuff, "Trp") == 0) {
			struct Amino *aminoTrp = (struct Amino*) malloc(sizeof(struct Amino));
			aminoTrp->name = "Trp";
			aminoTrp->trip = temps[incr].trip;
			aminoTrp->freq = temps[incr].freq;
			aminoTrp->rscu = 1;
			aminoTrp->count = 0;
			arrTrp[0] = aminoTrp;
		}
		else if(strcmp(tempBuff, "End") == 0) {
			struct Amino *aminoEnd = (struct Amino*) malloc(sizeof(struct Amino));
			aminoEnd->name = "End";
			aminoEnd->trip = temps[incr].trip;
			aminoEnd->freq = temps[incr].freq;
			aminoEnd->rscu = 0;
			aminoEnd->count = 0;
			arrEnd[endC] = aminoEnd;
			endC++;
		}
		else if(strcmp(tempBuff, "Cys") == 0) {
			struct Amino *aminoCys = (struct Amino*) malloc(sizeof(struct Amino));
			aminoCys->name = "Cys";
			aminoCys->trip = temps[incr].trip;
			aminoCys->freq = temps[incr].freq;
			aminoCys->rscu = 0;
			aminoCys->count = 0;
			arrCys[cysC] = aminoCys;
			cysC++;
		}
		else if(strcmp(tempBuff, "Tyr") == 0) {
			struct Amino *aminoTyr = (struct Amino*) malloc(sizeof(struct Amino));
			aminoTyr->name = "Tyr";
			aminoTyr->trip = temps[incr].trip;
			aminoTyr->freq = temps[incr].freq;
			aminoTyr->rscu = 0;
			aminoTyr->count = 0;
			arrTyr[tyrC] = aminoTyr;
			tyrC++;
		}
		else if(strcmp(tempBuff, "Leu") == 0) {
			struct Amino *aminoLeu = (struct Amino*) malloc(sizeof(struct Amino));
			aminoLeu->name = "Leu";
			aminoLeu->trip = temps[incr].trip;
			aminoLeu->freq = temps[incr].freq;
			aminoLeu->rscu = 0;
			aminoLeu->count = 0;
			arrLeu[leuC] = aminoLeu;
			leuC++;
		}
		else if(strcmp(tempBuff, "Phe") == 0) {
			struct Amino *aminoPhe = (struct Amino*) malloc(sizeof(struct Amino));
			aminoPhe->name = "Phe";
			aminoPhe->trip = temps[incr].trip;
			aminoPhe->freq = temps[incr].freq;
			aminoPhe->rscu = 0;
			aminoPhe->count = 0;
			arrPhe[pheC] = aminoPhe;
			pheC++;
		}
		else if(strcmp(tempBuff, "Gln") == 0) {
			struct Amino *aminoGln = (struct Amino*) malloc(sizeof(struct Amino));
			aminoGln->name = "Gln";
			aminoGln->trip = temps[incr].trip;
			aminoGln->freq = temps[incr].freq;
			aminoGln->rscu = 0;
			aminoGln->count = 0;
			arrGln[glnC] = aminoGln;
			glnC++;
		}
		else if(strcmp(tempBuff, "His") == 0) {
			struct Amino *aminoHis = (struct Amino*) malloc(sizeof(struct Amino));
			aminoHis->name = "His";
			aminoHis->trip = temps[incr].trip;
			aminoHis->freq = temps[incr].freq;
			aminoHis->rscu = 0;
			aminoHis->count = 0;
			arrHis[hisC] = aminoHis;
			hisC++;
		}
		else if(strcmp(tempBuff, "Pro") == 0) {
			struct Amino *aminoPro = (struct Amino*) malloc(sizeof(struct Amino));
			aminoPro->name = "Pro";
			aminoPro->trip = temps[incr].trip;
			aminoPro->freq = temps[incr].freq;
			aminoPro->rscu = 0;
			aminoPro->count = 0;
			arrPro[proC] = aminoPro;
			proC++;
		}
		else{
		}
		strcpy(tempBuff, rep);
	}

	strcpy(sequence, rep);
	strcpy(sequence, buffer);
	strcpy(other, rep);
	loop = 1;
	
	// count occurences of each individual codon.
	while(((strlen(sequence) - 1) % 3) == 0) {
		strncpy(other, sequence, 4);
		other[3] = '\0';

		if(strcmp(other, "") == 0) {
			break;
		}
		else {

			if(countCodon(arrGly, other, 4)) {
				strcpy(sequence, rep);
		                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
               			loop++;
		                strcpy(other, rep);
				continue;
			}
			else if(countCodon(arrGlu, other, 2)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrAsp, other, 2)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrVal, other, 4)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrAla, other, 4)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrArg, other, 6)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrSer, other, 6)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrLys, other, 2)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrAsn, other, 2)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrMet, other, 1)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrIle, other, 3)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrThr, other, 4)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrTrp, other, 1)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrEnd, other, 3)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrCys, other, 2)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrTyr, other, 2)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrLeu, other, 6)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrPhe, other, 2)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrGln, other, 2)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrHis, other, 2)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else if(countCodon(arrPro, other, 4)) {
                                strcpy(sequence, rep);
                                strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
                                loop++;
                                strcpy(other, rep);
                                continue;
                        }
			else {
			}
		}

		strcpy(sequence, rep);
		strncpy(sequence, &buffer[(3 * loop)], strlen(buffer));
		loop++;
		strcpy(other, rep);
	}

	// calculate rscu for each codon
	rscu(arrGly, 4);
	rscu(arrGlu, 2);
	rscu(arrAsp, 2);
	rscu(arrVal, 4);
	rscu(arrAla, 4);
	rscu(arrArg, 6);
	rscu(arrSer, 6);
	rscu(arrLys, 2);
	rscu(arrAsn, 2);
	rscu(arrMet, 1);
	rscu(arrIle, 3);
	rscu(arrThr, 4);
	rscu(arrTrp, 1);
	rscu(arrEnd, 3);
	rscu(arrCys, 2);
	rscu(arrTyr, 2);
	rscu(arrLeu, 6);
	rscu(arrPhe, 2);
	rscu(arrGln, 2);
	rscu(arrHis, 2);
	rscu(arrPro, 4);

	// calculate cai
	double caiVal = cai(arrGly, arrGlu, arrAsp, arrVal, arrAla, arrArg, arrSer, arrLys, arrAsn, arrMet, arrIle, arrThr, arrTrp, arrEnd, arrCys, arrTyr, arrLeu, arrPhe, arrGln, arrHis, arrPro, ((strlen(buffer) - 1) / 3), buffer);

	// the next 50 lines or so are the calculation for scuo
	double glyOi = calcOi(calcHi(arrGly, 4), 4);
	double gluOi = calcOi(calcHi(arrGlu, 2), 2);
	double aspOi = calcOi(calcHi(arrAsp, 2), 2);
	double valOi = calcOi(calcHi(arrVal, 4), 4);
	double alaOi = calcOi(calcHi(arrAla, 4), 4);
	double argOi = calcOi(calcHi(arrArg, 6), 6);
	double serOi = calcOi(calcHi(arrSer, 6), 6);
	double lysOi = calcOi(calcHi(arrLys, 2), 2);
	double asnOi = calcOi(calcHi(arrAsn, 2), 2);
	double ileOi = calcOi(calcHi(arrIle, 3), 3);
	double thrOi = calcOi(calcHi(arrThr, 4), 4);
	double endOi = calcOi(calcHi(arrEnd, 3), 3);
	double cysOi = calcOi(calcHi(arrCys, 2), 2);
	double tyrOi = calcOi(calcHi(arrTyr, 2), 2);
	double leuOi = calcOi(calcHi(arrLeu, 6), 6);
	double pheOi = calcOi(calcHi(arrPhe, 2), 2);
	double glnOi = calcOi(calcHi(arrGln, 2), 2);
	double hisOi = calcOi(calcHi(arrHis, 2), 2);
	double proOi = calcOi(calcHi(arrPro, 4), 4);

	double listOi[19] = { glyOi, gluOi, aspOi, valOi, alaOi, argOi, serOi, lysOi, asnOi, ileOi, thrOi, endOi, cysOi, tyrOi, leuOi, pheOi, glnOi, hisOi, proOi };

	double fiP1[19] = { fiDenom(arrGly, 4), fiDenom(arrGlu, 2), fiDenom(arrAsp, 2), fiDenom(arrVal, 4), fiDenom(arrAla, 4), fiDenom(arrArg, 6), fiDenom(arrSer, 6), fiDenom(arrLys, 2), fiDenom(arrAsn, 2), fiDenom(arrIle, 3), fiDenom(arrThr, 4), fiDenom(arrEnd, 3), fiDenom(arrCys, 2), fiDenom(arrTyr, 2), fiDenom(arrLeu, 6), fiDenom(arrPhe, 2), fiDenom(arrGln, 2), fiDenom(arrHis, 2), fiDenom(arrPro, 2) };

	double fi[19] = { 0 };

	int iii = 0;

	for(; iii < 19; iii++) {
		int jjj = 0;
		for(; jjj < 19; jjj++) {
			if(iii == jjj) {
				continue;
			}
			else {
				fi[iii] = fi[iii] + fiP1[jjj];
			}
		}
	}

	iii = 0;

	for(; iii < 19; iii++) {
		fi[iii] = fiP1[iii] / fi[iii];
	}

	iii = 0;
	double oFinal = 0;
	for(; iii < 19; iii++) {
		if(fi[iii] != 0) {
			fi[iii] = fi[iii] * listOi[iii];
		}
	}
	for(iii = 0; iii < 19; iii++) {
		oFinal = oFinal + fi[iii];
	}

	printf("Please choose which codon bias measure you would like to see. Type 'help' for a list of options, or type 'q' to exit the program.\n");

	// parse user input and display the requested information.
	while(true) {

		printf("\n > ");

		char input[20];

		scanf("%s", input);

		if(strcmp(input, "help") == 0 || strcmp(input, "Help") == 0 || strcmp(input, "HELP") == 0) {
			printf("\n   > 'rscu'\n\n     > RSCU stands for relative synonymous codon usage, and is a value between 0 and 6\n       (with the upper bound being the number of synonymous codons for the given amino\n       acid), which represents the over- or under-representation of a given codon\n       in a gene.\n\n   > 'cai'\n\n     > CAI stands for Codon Adaptation Index, and is a value between 0 and 1, which\n       represents the geometric mean of the weight (in codons) of each codon\n       over the sequence.\n\n   > 'scuo'\n\n     > SCUO stands for Synonymous Codon Usage Order and is a value between 0 and 1,\n       which represents the averaged normalized difference of the\n       codon distribution. The higher the value, the more biased on average are the\n       gene's codon selections.\n\n   > 'q' to quit\n");
			continue;
		}
		else if(strcmp(input, "rscu") == 0 || strcmp(input, "Rscu") == 0 || strcmp(input, "RSCU") == 0) {
			printf("\n");
			dispRscu(arrGly, arrGlu, arrAsp, arrVal, arrAla, arrArg, arrSer, arrLys, arrAsn, arrMet, arrIle, arrThr, arrTrp, arrEnd, arrCys, arrTyr, arrLeu, arrPhe, arrGln, arrHis, arrPro);
			continue;
		}
		else if(strcmp(input, "cai") == 0 || strcmp(input, "Cai") == 0 || strcmp(input, "CAI") == 0) {
			printf("\nThe CAI value is %f\n", caiVal);
		}
		else if(strcmp(input, "scuo") == 0 || strcmp(input, "Scuo") == 0 || strcmp(input, "SCUO") == 0) {
			printf("\nThe SCUO value is %f\n", oFinal);
		}
		else if(strcmp(input, "q") == 0 || strcmp(input, "Q") == 0) {
			printf("\n > Goodbye!\n\n");
			return 1;
		}
		else {
			printf("Unrecognized command. Please try again ('help' for a list of valid commands).");
			continue;
		}
	}
}

