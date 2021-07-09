//
// Created by Colin Small on 7/6/21.
//

#include <iostream>
#include <fstream>

#include "estimate_TFestimate.h"
#include "_helper.h"
#include "_parameter.h"
#include "_probability.h"
#include "_filter.h"

double MERGED_VAF_THRESHOLD;

int ZERO_MAPQUAL_COUNT_FOR_ESTIMATION;


void importPredefinedHotspotPrior(std::fstream file) {

    std::string line;
    std::vector<std::string> splitVector;
    while ( getline(file, line)) {
        splitVector = split(line, "\t");
        HOTSPOT_PRIOR.emplace_back(splitVector[0], std::stoi(splitVector[1]));
    }
}

void selectHotSpot(std::string chrom, std::string pos, std::string basestring, std::vector<double> quallist,
                   std::vector<double> maplist, std::string basestringNormal, std::vector<double> quallistNormal,
                   std::vector<double> maplistNormal, std::string basestringExtendedFrags,
                   std::vector<double> quallistExtendedFrags, std::vector<double> maplistExtendedFrags,
                   std::string basestringNotCombined, std::vector<double> quallistNotCombined,
                   std::vector<double> maplistNotCombined, std::map<char, int> basecountNotCombined,
                   std::map<char, int> basecountExtendedFrags, char variantBase, int nread) {

    if ( nread < DEPTH_FOR_ESTIMATION )
        return;

    std::string basestringNormalUpper = toUpper(basestringNormal);

    double basestringNormalVariantCount = std::count(basestringNormalUpper.begin(), basestringNormalUpper.end(),
                                                     variantBase);
    if ( basestringNormalVariantCount > GERMLINE_VARIANT_COUNT )
        return;

    std::map<char, int> basecount = countBase(basestring);
    std::string basestringMerge = basestringExtendedFrags + basestringNotCombined;

    std::vector<double> quallistMerge;
    quallistMerge.reserve(quallistExtendedFrags.size() + quallistNotCombined.size());
    quallistMerge.insert(quallistMerge.end(), quallistExtendedFrags.begin(), quallistNotCombined.end());
    quallistMerge.insert(quallistMerge.end(), quallistExtendedFrags.begin(), quallistNotCombined.end());

    std::vector<double> maplistMerge;
    maplistMerge.reserve(maplistExtendedFrags.size() + maplistNotCombined.size());
    maplistMerge.insert(maplistMerge.end(), maplistExtendedFrags.begin(), maplistExtendedFrags.end());
    maplistMerge.insert(maplistMerge.end(), maplistExtendedFrags.begin(), maplistExtendedFrags.end());

    // Require DEPTH_FOR_ESTIMATION fragments
    if ( basestringMerge.size() < DEPTH_FOR_ESTIMATION )
        return;
    // Exclude loci without enough germline information
    if ( basestringNormal.size() < DEPTH_FOR_GERMLINE )
        return;
    // Require high-quality variant alleles
    if ( basestring.size() == 0 || basestringMerge.size() == 0 || basestringNormal.size() == 0 )
        return;

    std::string basestringUpper = toUpper(basestring);
    int countNonRefRead = std::count(basestringUpper.begin(), basestringUpper.end(), variantBase);
    double meanQual = sumVector(quallist) / (double) quallist.size();
    if ( countNonRefRead < meanQual * (double) basestring.size())
        return;

    std::string baseStringMergeUpper = toUpper(basestringMerge);
    int countNonRefReadMerge = std::count(baseStringMergeUpper.begin(), baseStringMergeUpper.end(), variantBase);
    double meanQuallistMerge = sumVector(quallistMerge) / (double) quallistMerge.size();
    if ( countNonRefReadMerge < meanQuallistMerge * (double) basestringMerge.size())
        return;

    // Merged VAF lower than threshold
    // TODO: This corresponds to line 56 in p2.estimate.TFestimate.py, make sure that cnt_var_read_merge in that file
    // TODO: is the same as cnt_non_ref_read_merge and not a typo
    if ( countNonRefReadMerge / (double) basestringMerge.size() > MERGED_VAF_THRESHOLD )
        return;

    // Exclude possible germline variants
    if ((double) basestringNormalVariantCount / (double) basestringNormal.size() > GERMLINE_VAF_THRESHOLD_IN_NORMAL )
        return;

    // Require high maping quality in both tumor and normal
    int maplistZeroCount = std::count(maplist.begin(), maplist.end(), 0);
    int maplistNormalZeroCount = std::count(maplistNormal.begin(), maplistNormal.end(), 0);
    int maplistMergeZeroCount = std::count(maplistMerge.begin(), maplistMerge.end(), 0);
    if ( maplistZeroCount > ZERO_MAPQUAL_COUNT_FOR_ESTIMATION ||
         maplistNormalZeroCount > ZERO_MAPQUAL_COUNT_FOR_ESTIMATION ||
         maplistMergeZeroCount > ZERO_MAPQUAL_COUNT_FOR_ESTIMATION )
        return;

    std::vector<double> maplistGreaterThanZero = accessMultipleIndicesVector(maplist,
                                                                             getIndicesWhereGreaterThan(maplist, 0.0));
    std::vector<double> maplistNormalGreaterThanZero = accessMultipleIndicesVector(maplistNormal,
                                                                                   getIndicesWhereGreaterThan(
                                                                                           maplistNormal, 0.0));
    std::vector<double> maplistMergeGreaterThanZero = accessMultipleIndicesVector(maplistMerge,
                                                                                  getIndicesWhereGreaterThan(
                                                                                          maplistMerge, 0.0));
    double maplistMean = meanVector(maplistGreaterThanZero);
    double maplistNormalMean = meanVector(maplistNormalGreaterThanZero);
    double maplistMergeMean = meanVector(maplistMergeGreaterThanZero);
    if ( maplistMean > ZERO_MAPQUAL_COUNT_FOR_ESTIMATION ||
         maplistNormalMean > ZERO_MAPQUAL_COUNT_FOR_ESTIMATION ||
         maplistMergeMean > ZERO_MAPQUAL_COUNT_FOR_ESTIMATION )
        return;

    // Exclude strand bias
    char strandBiasMerge = filterStrandBiasMerge(basecountNotCombined, basecountExtendedFrags, variantBase);

    // Require enough percentage in non-reference reads
    std::tuple<char, char> bothObservedAboveAverageMerge = filterBothStrandAboveAverageMerge(basecountNotCombined,
                                                                                             basecountExtendedFrags);
    char bothObservedMerge = std::get<0>(bothObservedAboveAverageMerge);
    char aboveAverageMerge = std::get<1>(bothObservedAboveAverageMerge);

    char supportFragCountMerge = filterSupport(basestringMerge, variantBase);
    char strandBiasUnmerge = filterStrandBiasUnmerge(basecount, variantBase);

    // Require enough percentage in non-reference reads
    std::tuple<char, char> bothObservedAboveAverageUnmerge = filterBothStrandAboveAverageUnmerge(basecount);
    char bothObservedUnmerge = std::get<0>(bothObservedAboveAverageUnmerge);
    char aboveAverageUnmerge = std::get<1>(bothObservedAboveAverageUnmerge);

    char supportReadCountUnmerge = filterSupport(basestring, variantBase);

    if ( strandBiasMerge != 'T' or bothObservedMerge != 'T' or aboveAverageMerge != 'T' or
         supportFragCountMerge != 'T' )
        return;
    if ( strandBiasUnmerge != 'T' or bothObservedUnmerge != 'T' or aboveAverageUnmerge != 'T' or
         supportReadCountUnmerge != 'T' )
        return;

    HOTSPOT.emplace_back(hotspot(meanQuallistMerge, basestringNormalVariantCount,
                         (double) basestringNormalVariantCount /
                         basestringMerge.size(),
                         chrom, pos, basestring, quallist, maplist, basestringNormal, quallistNormal, maplistNormal,
                         basestringExtendedFrags, quallistExtendedFrags, maplistExtendedFrags, basestringNotCombined,
                         quallistNotCombined, maplistNotCombined, variantBase, nread));
}

// Unused in original python script
//void selectHotspotFromPredefinedHotspotPrior(std::string_view chrom, int pos, std::string_view basestring,
//                                             std::vector<double> & quallist, std::vector<double> & maplist,
//                                             std::string_view basestringNormal, std::vector<double> & quallistNormal,
//                                             std::vector<double> & maplistNormal,
//                                             std::string_view basestringExtendedFrags,
//                                             std::vector<double> & quallistExtendedFrags,
//                                             std::vector<double> & maplistExtendedFrags,
//                                             std::string_view basestringNotCombined,
//                                             std::vector<double> & quallistNotCombined,
//                                             std::vector<double> & maplistNotCombined,
//                                             std::map<char, int> & basecountNotCombined,
//                                             std::map<char, int> & basecountExtendedFrags, char variantBase) {
//    if( vecContains(HOTSPOT_PRIOR, std::tuple<std::string, int>(chrom, pos))){
//        std::cout << "(" << chrom << ", " << pos << ")";
//        selectHotSpot(chrom, pos, basestring, quallist, maplist, basestringNormal, quallistNormal, maplistNormal,  basestringExtendedFrags, quallistExtendedFrags, maplistExtendedFrags, basestringNotCombined, quallistNotCombined, maplistNotCombined, basecountNotCombined, basecountExtendedFrags, variantBase);
//    }
//}

double calculateTumorFractionLikelihood(double tumor_fraction, std::string basestring_merge,
                                        std::vector<double> & quallist_merge, std::vector<double> & maplist_merge,
                                        std::string basestring_normal, std::vector<double> & quallist_normal,
                                        std::vector<double> & maplist_normal, char variant_base) {

    double likelihood = 0.0;
    std::map<std::string, double>::iterator it;
    std::string index;
    double loglikelihood;
    for ( it = PRIOR.begin() ; it != PRIOR.end() ; it++ ) {
        index = it->first;
        loglikelihood =
                calculate_joint_genotype_tumor_fraction_loglikelihood(tumor_fraction, basestring_merge, quallist_merge,
                                                                      maplist_merge, index, variant_base) +
                calculate_joint_genotype_tumor_fraction_loglikelihood(0.0, basestring_normal, quallist_normal,
                                                                      maplist_normal, index, variant_base);
        likelihood += exp(loglikelihood + std::log(EST_PRIOR.find(index)->second));
    }
    return likelihood;
}

std::tuple<int, std::vector<double>> estimateTumorFraction(std::vector<hotspot> HOTSPOT) {

    std::vector<int> ratio;
    ratio.reserve(MAXSEARCH_WITH_NORMAL);
    for ( int i = 0 ; i < MAXSEARCH_WITH_NORMAL ; i++ ) {
        ratio[i] = GRIDWIDTH * i;
    }

    std::vector<double> tumorFractionLikelihood(0, MAXSEARCH_WITH_NORMAL);
    for ( int ratioind = 0 ; ratioind < MAXSEARCH_WITH_NORMAL ; ratioind++ ) {
        for ( hotspot spot : HOTSPOT ) {
            std::string basestringMerge = spot.basestringExtendedFrags + spot.basestringNotCombined;
            std::vector<double> quallistMerge = spot.quallistExtendedFrags;
            quallistMerge.insert(quallistMerge.begin(), spot.quallistNotCombined.begin(),
                                 spot.quallistNotCombined.end());
            std::vector<double> maplistMerge = spot.maplistExtendedFrags;
            maplistMerge.insert(maplistMerge.begin(), spot.maplistNotCombined.begin(), spot.maplistNotCombined.end());

            double likelihood = calculateTumorFractionLikelihood(ratio[ratioind], basestringMerge, quallistMerge,
                                                                 maplistMerge, spot.basestringNormal,
                                                                 spot.quallistNormal, spot.maplistNormal,
                                                                 spot.variantBase);
            tumorFractionLikelihood[ratioind] += std::log(likelihood) / (spot.basestringExtendedFrags.size() +
                                                                         spot.basestringNotCombined.size() +
                                                                         spot.basestringNormal.size());
        }
    }

    int estind = std::max_element(tumorFractionLikelihood.begin(), tumorFractionLikelihood.end()) -
                 tumorFractionLikelihood.begin();
    int est = ratio[estind];
    return std::tuple<int, std::vector<double>>(est, tumorFractionLikelihood);
}

int main(int argc, char *argv[]) {

    std::string filename = argv[1];
    MERGED_VAF_THRESHOLD = std::stod(argv[2]);
    std::string file_prefix = argv[3];
    double depth = std::stod(argv[4]);
    std::string VAF_output = argv[5];
    std::string estimate_output = argv[6];
    GERMLINE_VARIANT_COUNT = getGermlineVariantCountThreshold(depth);

    std::ifstream file(filename);

    std::string lineText;
    while ( std::getline(file, lineText)) {
        std::vector<std::string> splitString = split(lineText, "\t");
        std::string chrom = splitString[0];

        if ( chrom.find("X") == std::string::npos ||
             chrom.find("M") == std::string::npos ||
             chrom.find("Y") == std::string::npos )
            continue;

        std::string pos = splitString[1];
        int nread = splitString[19].size();
        std::string basestringAll = splitString[4];
        std::string basestring = splitString[19];

        if ( basestring.size() == 0 )
            continue;

        std::vector<double> quallist = stringToQual(splitString[20]);
        std::vector<double> maplist = stringToQual(splitString[21]);

        if ( splitString.size() < 31 )
            continue;

        std::string basestring_normal_all = splitString[8];
        std::string basestring_normal = splitString[22];
        std::vector<double> quallist_normal = stringToQual(splitString[23]);
        std::vector<double> maplist_normal = stringToQual(splitString[24]);
        std::string basestring_extendedFrags = splitString[25];
        std::vector<double> quallist_extendedFrags = stringToQual(splitString[26]);
        std::vector<double> maplist_extendedFrags = stringToQual(splitString[27]);
        std::string basestring_notCombined = splitString[28];
        std::vector<double> quallist_notCombined = stringToQual(splitString[29]);
        std::vector<double> maplist_notCombined = stringToQual(splitString[30]);
        std::map<char, int> basecount_notCombined_all = countBase(splitString[12]);
        std::map<char, int> basecount_extendedFrags_all = countBase(splitString[16]);
        char variantBaseAll = findMajorVariant(basecount_notCombined_all, basecount_extendedFrags_all);
        std::map<char, int> basecountNotCombined = countBase(basestring_notCombined);
        std::map<char, int> basecountExtendedFrags = countBase(basestring_extendedFrags);
        char variantBase = findMajorVariant(basecountNotCombined, basecountExtendedFrags);
        std::map<char, int> basecountTumorAll = countBase(basestringAll);
        std::map<char, int> basecountTumor = countBase(basestring);
        char trialleleTumorAll = filterTriallelicPosition(basecountTumorAll, variantBaseAll, depth);
        char trialleleTumor = filterTriallelicPosition(basecountTumor, variantBaseAll, depth);

        if ( variantBase != variantBaseAll )
            continue;
        if ( std::count(basestring.begin(), basestring.end(), std::toupper(variantBase)) +
             std::count(basestring.begin(), basestring.end(), std::tolower(variantBase)) == 0 )
            continue;
        if ( trialleleTumorAll == 'F' or trialleleTumor == 'F' )
            continue;

        std::string basestringNormallAllUpper = toUpper(basestring_normal_all);
        if ( std::count(basestringNormallAllUpper.begin(), basestringNormallAllUpper.end(), variantBase) >
             GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER )
            continue;
        if ( std::count(basestringNormallAllUpper.begin(), basestringNormallAllUpper.end(), variantBaseAll) >
             GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER )
            continue;
        selectHotSpot(chrom, pos, basestring, quallist, maplist, basestring_normal, quallist_normal, maplist_normal,
                      basestring_extendedFrags, quallist_extendedFrags, maplist_extendedFrags, basestring_notCombined,
                      quallist_notCombined, maplist_notCombined, basecountNotCombined, basecountExtendedFrags,
                      variantBase, nread);
    }
    file.close();

    if ( HOTSPOT.size() > MAX_NUMBER_OF_HOTPOST ) {
//        from operator import itemgetter
//
//                a = sorted(HOTSPOT, key=itemgetter(0))
//        a = sorted(a, key=itemgetter(2), reverse = True)
//        a = sorted(a, key=itemgetter(1))
//        HOTSPOT = a[0:MAX_NUMBER_OF_HOTSPOT]
    }

    std::tuple<int, std::vector<double>> estTumorFractionLikelihood = estimateTumorFraction(HOTSPOT);
    int est = std::get<0>(estTumorFractionLikelihood);
    std::vector<double> tumorFractionLikelihood = std::get<1>(estTumorFractionLikelihood);

    std::vector<double> VAF;
    VAF.reserve(HOTSPOT.size());
    for ( hotspot & hs : HOTSPOT ) {
        std::string basestringCombined = hs.basestringNotCombined + hs.basestringExtendedFrags;
        int basecount = std::count(basestringCombined.begin(), basestringCombined.end(), std::toupper(hs.variantBase)) +
                        std::count(basestringCombined.begin(), basestringCombined.end(), std::tolower(hs.variantBase));
        double vaf = (double) basecount / basestringCombined.size();
        VAF.push_back(vaf);
    }
    std::sort(VAF.begin(), VAF.end());

    std::ofstream vafOutputFile(VAF_output);
    for ( double & vaf : VAF )
        vafOutputFile << vaf << '\n';
    vafOutputFile.close();

    std::ofstream estimateOutputFile(estimate_output);
    estimateOutputFile << est << '\n';
    estimateOutputFile << "====================" << '\n';
    for ( double & tfl : tumorFractionLikelihood )
        estimateOutputFile << tfl << '\n';
    estimateOutputFile.close();
    return 0;
}