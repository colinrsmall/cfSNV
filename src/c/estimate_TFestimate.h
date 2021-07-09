//
// Created by Colin Small on 7/6/21.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "hotspot.h"

#ifndef CFSNV_ESTIMATE_TFESTIMATE_H
#define CFSNV_ESTIMATE_TFESTIMATE_H

void importPredefinedHotspotPrior(std::fstream file);

void selectHotSpot(std::string chrom,
                   std::string pos,
                   std::string basestring,
                   std::vector<double> quallist,
                   std::vector<double> maplist,
                   std::string basestringNormal,
                   std::vector<double> quallistNormal,
                   std::vector<double> maplistNormal,
                   std::string basestringExtendedFrags,
                   std::vector<double> quallistExtendedFrags,
                   std::vector<double> maplistExtendedFrags,
                   std::string basestringNotCombined,
                   std::vector<double> quallistNotCombined,
                   std::vector<double> maplistNotCombined,
                   std::map<char, int> basecountNotCombined,
                   std::map<char, int> basecountExtendedFrags,
                   char variantBase, int nread);

void selectHotspotFromPredefinedHotspotPrior(std::string_view chrom,
                                             std::string_view pos,
                                             std::string_view basestring,
                                             std::vector<double>& quallist,
                                             std::vector<double>& maplist,
                                             std::string_view basestringNormal,
                                             std::vector<double>& quallistNormal,
                                             std::vector<double>& maplistNormal,
                                             std::string_view basestringExtendedFrags,
                                             std::vector<double>& quallistExtendedFrags,
                                             std::vector<double>& maplistExtendedFrags,
                                             std::string_view basestringNotCombined,
                                             std::vector<double>& quallistNotCombined,
                                             std::vector<double>& maplistNotCombined,
                                             std::map<char, int>& basecountNotCombined,
                                             std::map<char, int>& basecountExtendedFrags,
                                             char variantBase);

int calculateTumorFractionLikelihood(double tumor_fraction,
                                     std::string_view basestring_merge,
                                     std::vector<double>& quallist_merge,
                                     std::vector<double>& maplist_merge,
                                     std::string_view basestring_normal,
                                     std::vector<double>& quallist_normal,
                                     std::vector<double>& maplist_normal,
                                     char variant_base);

std::tuple<int, std::vector<double>> estimateTumorFraction(hotspot hs);

#endif //CFSNV_ESTIMATE_TFESTIMATE_H
