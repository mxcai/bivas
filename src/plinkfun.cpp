//
//  plinkfun.cpp
//  bivas_Xcode
//
//  Created by CAI Mingxuan on 2017/3/28.
//  Copyright © 2017年 CAI Mingxuan. All rights reserved.
//

#include "plinkfun.hpp"
#include <stdio.h>
#include <math.h>

int getLineNum(string filename){
    FILE *pf = fopen(filename.c_str(), "r");
    char buf[10000];
    int lineCnt = 0;
    if (!pf)
        return -1;
    while (fgets(buf, 10000, pf))
        lineCnt++;
    fclose(pf);
    return lineCnt;
}

void getFourGentype(int* geno, std::bitset<8> bits){
    int idx = 0;
    for (int j=0; j < 8; j = j + 2) {
        if(bits[j] && bits[j+1]){
            geno[idx] = 0;
        }else if(!bits[j] && !bits[j+1]){
            geno[idx] = 2;
        }else if(!bits[j] && bits[j+1]){
            geno[idx] = 1;
        }else if(bits[j] && !bits[j+1]){
            geno[idx] = 3;
        }
        idx++;
    }
}


void readPlink(string stringname,int N, int P, int* X){

    FILE *fp;
    unsigned char buff[3];
    string bedfile = stringname +".bed";
    fp = fopen(bedfile.c_str(), "rb");
    if (!fp) return;
    fread(buff, sizeof(char), 3, fp);

    std::bitset<8> magic1(buff[0]);
    std::bitset<8> magic2(buff[1]);
    std::bitset<8> mode0(buff[2]);

    if(magic1.to_ulong() != 108 || magic2.to_ulong() != 27){
    }

    unsigned long mode =  mode0.to_ulong();

    int n = 0;
    long charNum = ceil(N*1.0/4)*10000;
    int leftGenoNum = ceil(N*1.0/4)*P;
    int nblock = ceil(N*1.0/4);
    int nSNP = 0;
    while (!feof(fp)) {
        if(leftGenoNum <= 0)
            break;
        if(leftGenoNum <= charNum){
            charNum  = leftGenoNum;
        }
        char* genotype = new char[charNum];
        fread(genotype, sizeof(char), charNum, fp);
        int* geno = new int[4];
        int nSNPc = int(charNum / nblock); //number of SNPs of this iteration
        int idx = 0;
        for (int i=0; i < nSNPc; i++) {
            for(int j=0; j < nblock - 1; j++){
                std::bitset<8> bits(genotype[idx]);
                getFourGentype(geno,bits);
                memcpy(X + nSNP * N + j*4, geno, 4*sizeof(int));
                idx++;
                leftGenoNum -= 1;
            }
            int left = N - (nblock - 1)*4;
            std::bitset<8> bits(genotype[idx]);
            getFourGentype(geno,bits);
            memcpy(X + nSNP * N + (nblock - 1)*4, geno, left*sizeof(int));
            idx++;
            leftGenoNum -= 1;
            nSNP ++;
        }
        delete[] geno;
        n++;
    }


}
