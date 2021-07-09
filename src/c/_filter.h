//
// Created by Colin Small on 7/7/21.
//

#ifndef CFSNV__FILTER_H
#define CFSNV__FILTER_H

#include "string"
#include "map"

char filterStrandBiasMerge(std::map<char, int> & basecountNotCombined,
                           std::map<char, int> & basecountExtendedFrags,
                           char variantBase);

char filterStrandBiasUnmerge(std::map<char, int> & basecount,
                             char variantBase);

std::tuple<char, char> filterBothStrandAboveAverageMerge(std::map<char, int> & basecountNotCombined,
                                       std::map<char, int> & basecountExtendedFrags);


char filterSupport(std::string basestringMerge, char variantBase);

std::tuple<char, char> filterBothStrandAboveAverageUnmerge(std::map<char, int> basecount);

int getGermlineVariantCountThreshold(double depth);

char findMajorVariant(std::map<char, int> & basecountNC, std::map<char, int> & basecountOF);

char filterTriallelicPosition(std::map<char, int> basecount, char var, int depth);
#endif //CFSNV__FILTER_H
