import sys
target_bed = open(sys.argv[1])
base_mapped = open(sys.argv[2])
n_read_base = float(base_mapped.readline().strip())
n_region_base = 0
for line in target_bed:
	sp = line.strip().split('\t')
	n_region_base = n_region_base + int(sp[2]) - int(sp[1])

base_mapped.close()
target_bed.close()
output = open(sys.argv[3], 'w')
output.write(str(float(n_read_base)/float(n_region_base) * 0.8) + '\n')
output.close()

