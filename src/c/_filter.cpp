//
// Created by Colin Small on 7/7/21.
//

#include "_filter.h"
#include "random"
#include "_parameter.h"
#include "iostream"
#include "_probability.h"

char filterStrandBiasMerge(std::map<char, int> & basecountNotCombined,
                                  std::map<char, int> & basecountExtendedFrags,
                                  char variantBase){
    int forward = basecountNotCombined['R'] + basecountNotCombined['A'] + basecountNotCombined['C'] + basecountNotCombined['T'] + basecountNotCombined['G'] + (basecountExtendedFrags['r'] + basecountExtendedFrags['a'] + basecountExtendedFrags['c'] + basecountExtendedFrags['t'] + basecountExtendedFrags['g'] + basecountExtendedFrags['R'] + basecountExtendedFrags['A'] + basecountExtendedFrags['C'] + basecountExtendedFrags['T'] + basecountExtendedFrags['G']);
    int reverse = basecountNotCombined['r'] + basecountNotCombined['a'] + basecountNotCombined['c'] + basecountNotCombined['t'] + basecountNotCombined['g'] + (basecountExtendedFrags['R'] + basecountExtendedFrags['A'] + basecountExtendedFrags['C'] + basecountExtendedFrags['T'] + basecountExtendedFrags['G'] + basecountExtendedFrags['r'] + basecountExtendedFrags['a'] + basecountExtendedFrags['c'] + basecountExtendedFrags['t'] + basecountExtendedFrags['g']);
    int forwardVar = basecountNotCombined[variantBase] + basecountExtendedFrags[std::tolower(variantBase)] + basecountExtendedFrags[variantBase];
    int reverseVar = basecountNotCombined[std::tolower(variantBase)] + basecountExtendedFrags[variantBase] + basecountExtendedFrags[std::tolower(variantBase)];

    if(reverseVar == 0 && forwardVar == 0)
        return 'F';
    if(reverseVar == 0 && reverse == 0)
        return 'T';
    if(reverseVar == 0 and reverse != 0) {
        // TODO: Make sure this actually calculates the same thing as the python script
        if ( binomialPMF(forwardVar, forwardVar + reverseVar, (double) reverse / (forward + reverse)) <
             STRAND_BIAS_BINOMIAL_PROB )
            return 'F';
        else
            return 'T';
    }
    else{
        double varRatio = (double)forwardVar/reverseVar;
        double ratio = (double)forward/reverse;
        if(varRatio/ratio > STRAND_BIAS_RATIO_VARIANT_TO_ALL ||
           ratio/varRatio > STRAND_BIAS_RATIO_VARIANT_TO_ALL)
            return 'F';
        else
            return 'T';
    }
}

std::tuple<char, char> filterBothStrandAboveAverageMerge(std::map<char, int> & basecountNotCombined,
                                       std::map<char, int> & basecountExtendedFrags)
{
    char bothObserved = 'F';
    char aboveAverage = 'F';

    for( char base : {'A', 'C', 'T', 'G'})
    {
        if((basecountNotCombined[std::tolower(base)] > 0 && basecountNotCombined[base] > 0 ) ||
           basecountExtendedFrags[std::tolower(base)] > 0 ||
           basecountExtendedFrags[base] > 0)
            bothObserved = 'T';

        if( basecountNotCombined[base] + basecountNotCombined[std::tolower(base)] + 2*basecountExtendedFrags[std::tolower(base)] + 2*basecountExtendedFrags[base] >=
            THRESHOLD_VARIANT_ALLELE_PROPORTION_IN_NONREF * (sumMap(basecountNotCombined) + 2*sumMap(basecountExtendedFrags) - basecountNotCombined['R'] - basecountNotCombined['r'] - 2*basecountExtendedFrags['R'] - 2*basecountExtendedFrags['r']))
            aboveAverage = 'T';
    }

    return std::tuple<char, char>(bothObserved, aboveAverage);
}

char filterSupport(std::string basestring, char variantBase){
    int variantCount = std::count(basestring.begin(), basestring.end(), variantBase);
    int variantCountLower = std::count(basestring.begin(), basestring.end(), std::tolower(variantBase));
    if( variantCount + variantCountLower < COUNT_FOR_STRONG_EVIDENCE )
        return 'F';
    else
        return 'T';
}

char filterStrandBiasUnmerge(std::map<char, int> & basecount,
                             char variantBase)
{
    int forward = basecount['R'] + basecount['A'] + basecount['C'] + basecount['T'] + basecount['G'];
    int reverse = basecount['r'] + basecount['a'] + basecount['c'] + basecount['t'] + basecount['g'];
    int forwardVar = basecount[variantBase];
    int reverseVar = basecount[std::tolower(variantBase)];

    if(reverseVar == 0 and reverse == 0)
        return 'T';
    if(reverseVar == 0 and reverse != 0)
    {
        if(binomialPMF(forwardVar, forwardVar+reverseVar, (double)reverse/(forward+reverse)) < STRAND_BIAS_BINOMIAL_PROB)
            return 'F';
        else
            return 'T';
    }
    if(forwardVar == 0 and forward == 0)
        return 'T';
    if(forwardVar == 0 and forward != 0)
    {
        if( binomialPMF(reverseVar, forwardVar+reverseVar, (double)reverse/(forward+reverse)) < STRAND_BIAS_BINOMIAL_PROB)
            return 'F';
        else
            return 'T';
    }
    else{
        double varRatio = (double)forwardVar/reverseVar;
        double ratio = (double)forward/reverse;
        if(varRatio/ratio > STRAND_BIAS_RATIO_VARIANT_TO_ALL ||
           ratio/varRatio > STRAND_BIAS_RATIO_VARIANT_TO_ALL)
            return 'F';
        else
            return 'T';
    }
}

std::tuple<char, char> filterBothStrandAboveAverageUnmerge(std::map<char, int> basecount){
    char bothObserved = 'F';
    char aboveAverage = 'T';

    for(char base : {'A', 'C', 'T', 'G'})
    {
        if(basecount[std::tolower(base)] > 0 and basecount[base] > 0)
            bothObserved = 'T';
        if(basecount[base] + basecount[std::tolower(base)] >= THRESHOLD_VARIANT_ALLELE_PROPORTION_IN_NONREF * (sumMap(basecount) - basecount['R'] - basecount['r']))
            aboveAverage = 'T';
    }

    return std::tuple<char, char>(bothObserved, aboveAverage);
}

int getGermlineVariantCountThreshold(double depth){
    return (int)std::log10(depth)-1;
}

char findMajorVariant(std::map<char, int> & basecountNC, std::map<char, int> & basecountOF){
    std::vector<char> alternativeNucleotide{'A', 'C', 'T', 'G'};
    std::vector<int> count;
    count.reserve(alternativeNucleotide.size());
    for(char& nucleotide : alternativeNucleotide){
        count.push_back(basecountNC[nucleotide] + basecountNC[std::tolower(nucleotide)] + 2*basecountOF[nucleotide] + 2*basecountOF[std::tolower(nucleotide)]);
    }
    int variantBaseIndex = std::max_element(count.begin(), count.end()) - count.begin();
    return alternativeNucleotide[variantBaseIndex];
}

char filterTriallelicPosition(std::map<char, int> basecount, char var, int depth){
    int varCount = basecount[var] + basecount[std::tolower(var)];
    int totalCount = sumMap(basecount);
//    if(totalCount == 0)
//        std::cout << totalCount << "," << varCount << "," << basecount << ',' << var;
    double varVAF = (double)varCount/totalCount;
    TRIALLELE_COUNT = (int)std::log10(depth) + 2;
    int otherCount = 0;

    for(char c : "ACGT"){
        if(c == var)
            continue;
        int tmpCount = basecount[c] + basecount[std::tolower(c)];
        if( tmpCount > otherCount)
            otherCount = tmpCount;
    }

    double otherVAF = (double)otherCount/totalCount;
    if(varCount == 0 or otherVAF/varVAF > TRIALLELE_VAF_RATIO)
        return 'F';
    if(otherVAF > TRIALLELE_VAF and totalCount >= 0)
        return 'F';
    if(otherCount > TRIALLELE_COUNT and totalCount <= 140)
        return 'F';
    return 'T';
}

double getNormalCountBinomThreshold(double depth){
    if(depth < 1000)
        return 0.05;
    else
        return 0.2;
}