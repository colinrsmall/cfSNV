from parameter import *
import numpy as np
import time

def string_to_base(basestring):
# input: a base string
# output: list of bases in the base string
# string type in python is iterable, so output is the same as input here
	return basestring

def string_to_qual(qualitystring):
# input: a qual string
# output: a list of error probabilities
# qual is the phred score of error probability
	prob = np.array([10.0**(-double(ord(i)-33)/10.0) for i in qualitystring])
	return prob


def filter_low_qual(basestring, quallist, maplist):
	id = ( quallist < BASEQUAL_THRESHOLD ) * ( maplist < MAPQUAL_THRESHOLD )
	new_basestring = "".join([basestring[i] for i in range(len(basestring)) if id[i] == True])
	new_quallist = quallist[id]
	new_maplist = maplist[id]
	return new_basestring, new_quallist, new_maplist

def filter_low_qual_with_string(basestring, quallist, maplist, qualstring, mapstring):
	id = ( quallist < BASEQUAL_THRESHOLD ) * ( maplist < MAPQUAL_THRESHOLD )
	new_basestring = "".join([basestring[i] for i in range(len(basestring)) if id[i] == True])
	new_quallist = quallist[id]
	new_maplist = maplist[id]
	new_qualstring = "".join([qualstring[i] for i in range(len(qualstring)) if id[i] == True])
	new_mapstring = "".join([mapstring[i] for i in range(len(mapstring)) if id[i] == True])
	return new_basestring, new_quallist, new_maplist, new_qualstring, new_mapstring



def find_major_variant(basecount_nc, basecount_ef):
# input: converted basestring
# output: variant with the max count of observations
#       haven't thought of the case where the max count corresponds to a variant with strong strand bias
	alternative_nucleotide = ['A', 'C', 'T', 'G']
	count = []
	for nucleotide in alternative_nucleotide:
		count.append( basecount_nc[nucleotide] + basecount_nc[nucleotide.lower()] + 2*basecount_ef[nucleotide] + 2*basecount_ef[nucleotide.lower()] )
	max_count = max(count)
	variant_base = alternative_nucleotide[count.index(max_count)]
	return variant_base


def observe_variant_probability(tumor_fraction, joint_genotype):
# input: tumor fraction, joint genotype
# output: theoretical allele frequencies of reference allele and non-reference allele
# combine tumor fraction and joint genotype to the frequencies of nucleotides
	variant_allele_observed_probability = dict([])
	if joint_genotype == 'AA/AA':
		variant_allele_observed_probability['R'] = 1.0
		variant_allele_observed_probability['X'] = 0.0
	elif joint_genotype == 'AA/AB':
		variant_allele_observed_probability['R'] = 1.0 - 0.5*tumor_fraction
		variant_allele_observed_probability['X'] = 0.5*tumor_fraction
	elif joint_genotype == 'AA/BB':
		variant_allele_observed_probability['R'] = 1.0 - tumor_fraction
		variant_allele_observed_probability['X'] = tumor_fraction
	elif joint_genotype == 'AB/AA':
		variant_allele_observed_probability['R'] = 0.5 + 0.5*tumor_fraction
		variant_allele_observed_probability['X'] = 0.5 - 0.5*tumor_fraction
	elif joint_genotype == 'AB/AB':
		variant_allele_observed_probability['R'] = 0.5
		variant_allele_observed_probability['X'] = 0.5
	elif joint_genotype == 'AB/BB':
		variant_allele_observed_probability['R'] = 0.5 - 0.5*tumor_fraction
		variant_allele_observed_probability['X'] = 0.5 + 0.5*tumor_fraction
	elif joint_genotype == 'BB/AA':
		variant_allele_observed_probability['R'] = tumor_fraction
		variant_allele_observed_probability['X'] = 1.0 - tumor_fraction
	elif joint_genotype == 'BB/AB':
		variant_allele_observed_probability['R'] = 0.5*tumor_fraction
		variant_allele_observed_probability['X'] = 1.0 - 0.5*tumor_fraction
	elif joint_genotype == 'BB/BB':
		variant_allele_observed_probability['R'] = 0.0
		variant_allele_observed_probability['X'] = 1.0
	else:
		print "error with input genotype"
	return variant_allele_observed_probability


def translate_basestring(basestring, variant_base):
# input: basestring in rebase file
# output: non-reference alleles, converted basestring, base count string
# find all observed non-reference alleles
# count occurrence of non-reference alleles
# convert reference alleles to 'R', and non-reference alleles to 'X'
	newstring = basestring
	newstring = newstring.replace(variant_base, 'X')
	newstring = newstring.replace(variant_base.lower(), 'X')
	newstring = newstring.replace('r', 'R')
	return np.array(list(newstring), dtype = str)


def count_base(basestring):
# input: basestring in rebase file
# output: base count dict
# count occurrence of non-reference alleles
	basecount = {}
	for i in "AaCcGgTtRr":
		basecount[i] = basestring.count(i)
	return basecount



def convert_to_Decimal(x):
# input: integer or double
# output: decimal object
	return Decimal(str(x))



def calculate_joint_genotype_tumor_fraction_loglikelihood(tumor_fraction, basestring, quallist, maplist, joint_genotype, variant_base):
# input: tumor fraction, converted basestring, base quality list, joint genotype
# output: log likelihood log( P( X | theta, joint_genotype ) )
	newbasestring = translate_basestring(basestring, variant_base)
	interest_id = (newbasestring == 'R') + (newbasestring == 'X')
	VAF = observe_variant_probability(tumor_fraction, joint_genotype)
	#print VAF, newbasestring, basestring, variant_base,  interest_id
#	print len(basestring), len(newbasestring), len(quallist), len(maplist)
	interest_quallist = quallist[interest_id]
	interest_maplist = maplist[interest_id]
	interest_VAFlist = np.array([VAF[i] for i in newbasestring[interest_id]])
	#print interest_VAFlist, interest_quallist, interest_maplist
	loglist = np.log(interest_VAFlist*(1-interest_quallist)*(1-interest_maplist) + (1-interest_VAFlist)*(interest_quallist+interest_maplist+interest_quallist*interest_maplist))
	loglikelihood = sum(loglist)
	return loglikelihood

def main():
	start = time.time()
	basestring = "rrRrrRRRrRRrrRRRRrrrrRrRrrRrRRrRaaaAAAAAAaaAaAaaaAaARrRRRrrrRrRrrRrRRRRrRrAAaaaaAaaaAAaaAaa"
	joint_genome = "BB/BB"
	maplist = np.array([0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.251188644, 0.000001, 0.000001, 0.000001])
	quallist = np.array([0.022535235, 0.00528222, 0.009541663, 0.013264347, 0.011090627, 0.003555052, 0.027080566, 0.015750399, 0.022683902, 0.031131367, 0.018204997, 0.00997583, 0.005136876, 0.019508333, 0.028207442, 0.028061324, 0.027953896, 0.028271181, 0.013522274, 0.026801979, 0.0016488, 0.005067841, 0.031046452, 0.019417774, 0.027828293, 0.023785416, 0.019658524, 0.006075595, 0.007924912, 0.029089328, 0.010350578, 0.017759347, 0.007546654, 0.021277069, 0.018700052, 0.006809801, 0.027671635, 0.006925604, 0.028333594, 0.024216763, 0.028855102, 0.029491, 0.011219381, 0.022957093, 0.014945312, 0.011761302, 0.001022148, 0.019425203, 0.026984942, 0.005588313, 0.028931969, 0.014299779, 0.00322236, 0.021123591, 0.004383596, 0.02223039, 0.025926345, 0.020896607, 0.019617726, 0.005057836, 0.012143053, 0.026517415, 0.020205884, 0.00856174, 0.007850494, 0.019941156, 0.028486627, 0.020179196, 0.013996312, 0.026127624, 0.01342404, 0.012679071, 0.011769746, 0.021126596, 0.013215724, 0.00816966, 0.018582238, 0.004288396, 0.014510591, 0.022444977, 0.020117276, 0.005726826, 0.003265246, 0.023436632, 0.024888499, 0.028886556, 0.026453318, 0.023166153, 0.017933538, 0.024731459, 0.03156874])
	tumor_fraction = 0.0
	variant_base = "A"
	for i in range(10000):
		calculate_joint_genotype_tumor_fraction_loglikelihood(tumor_fraction, basestring, quallist, maplist, joint_genome, variant_base)
	done = time.time()
	print done - start

main()