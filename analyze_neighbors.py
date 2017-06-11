import csv
import numpy as np
import random as r

linc_hits_file = 'top_hits.csv'
lincs_file = 'combined_lincRNAs_final.csv'
genes_file = 'unique_isoforms_NM_hg19.csv'
sam_hits_file = 'SAM_hits.csv'
baseline_file = '161218_baseline.csv'
sam_output_file = 'sam_output.csv'
rnaseq_output_file = 'rnaseq_output.csv'

power_output_file = 'power_output.csv'

bp_cutoff = 1000000
rnaseq_lincs = ['TCONS_00011252', 'NR_034078', 'TCONS_00010506', 'TCONS_00026344', 'TCONS_00015940', 'NR_125939', 'TCONS_00009861', 'NR_109890',
	'NR_033834', 'TCONS_00006579', 'NR_026873']




with open(lincs_file, 'rb') as gf:
	data = [row for row in csv.reader(gf.read().splitlines())]
	lincs = dict([(row[0], row[1:]) for i,row in enumerate(data)])

with open(genes_file, 'rb') as gf:
	data = [row for row in csv.reader(gf.read().splitlines())]
	genes = dict([(row[0], row[1:5]) for i,row in enumerate(data)])
	genes_map = dict([(row[0], row[5]) for i,row in enumerate(data)])

with open(linc_hits_file, 'rb') as gf:
	linc_hits = [row[0] for row in csv.reader(gf.read().splitlines())]

with open(sam_hits_file, 'rb') as gf:
	data = [row for row in csv.reader(gf.read().splitlines())]
	sam_hits = dict([(row[0], float(row[1])) for i,row in enumerate(data)])

# with open(baseline_file, 'rb') as gf:
# 	data = [row for row in csv.reader(gf.read().splitlines())]
# 	baseline = dict([(row[0], row[1:]) for i,row in enumerate(data[1:])])

def get_TSS(loc_list):
	strand = loc_list[1]
	if strand == "+":
		return loc_list[0], int(loc_list[2])
	elif strand == "-":
		return loc_list[0], int(loc_list[3])

def get_activation(exp, ctrl):
	exp = [float(e) for e in exp]
	ctrl = [float(c) for c in ctrl]
	return np.mean(exp)/(np.mean(ctrl)+0.0000001) > 1

with open(power_output_file, 'wb') as csvfile:
	csvwriter = csv.writer(csvfile)

	for l in ['TCONS_00006579']:
		power_file = l + '_power.csv'
		with open(power_file, 'rb') as gf:
			data = [row for row in csv.reader(gf.read().splitlines())]
			power_genes = dict([(row[0], row[1:]) for i,row in enumerate(data[1:])])
			print len(power_genes.keys())
		linc_chr, linc_loc = get_TSS(lincs[l])
		for g in genes.keys():
			gene_chr, gene_loc = get_TSS(genes[g])

			if gene_chr == linc_chr:
				if abs(linc_loc - gene_loc) < bp_cutoff:
					g = genes_map[g]
					if g in power_genes.keys():
						csvwriter.writerow([g] + power_genes[g])
					else:
						print g
		# 			# if g in sam_hits.keys():
		# 			# 	sam_csvwriter.writerow([l, g, genes_map[g], sam_hits[g]])
		# 			if g in diff_genes.keys():
		# 				# rnaseq_csvwriter.writerow([l, g] + diff_genes[g])
		# 				if l == 'TCONS_00026344':
		# 					a = get_activation(diff_genes[g][:7], diff_genes[g][7:])
		# 				else:
		# 					a = get_activation(diff_genes[g][:9], diff_genes[g][9:])
		# 				# if a:
		# 				if g not in sig_genes: sig_genes.append(g)
		# 			else:
		# 				if g == 'C6orf195': g = 'LINC01600'
		# 				if g == 'FLJ44313': continue
		# 				if g == 'LOC101927322': g = 'NM_001289967'
		# 				if g == 'CCDC11': g = 'CFAP53'
		# 				if g == 'KCNA6': continue
		# 				if g == 'GPR133': g = 'ADGRD1'
		# 				num_expressed = [0,0,0,0,0]
		# 				total_exp = 0
		# 				for i, exp in enumerate(baseline[g]):
		# 					total_exp += float(exp)
		# 					if float(exp) >= 1:
		# 						num_expressed[i/3] += 1
		# 				if (num_expressed[3] + num_expressed[4]) == 0 and (num_expressed[0] == 0 or num_expressed[1] == 0 or num_expressed[2] == 0):
		# 					if g not in undetected_genes and g not in lowly_exp_genes:
		# 						if total_exp < 0.001:
		# 							undetected_genes.append(g)
		# 						else:
		# 							lowly_exp_genes.append(g)
		# 					# rnaseq_csvwriter.writerow([l, g] + baseline[g] + ['undetected'])
		# 			if g not in all_genes: all_genes.append(g)
		# N = 10000
		# prob = float(n_diff)/len(genes_map.values())
		# count = 0
		# # if l == 'TCONS_00015940': sig_genes = sig_genes + ['MOB3B', 'IFNK', 'EQTN']
		# for i in range(N):
		# 	k = 0
		# 	for j in range(len(all_genes)):
		# 		if r.random() <= prob: k += 1
		# 	if k == len(sig_genes):
		# 		count += 1

		# # if l in rnaseq_lincs:
		# print l, len(diff_genes.keys()), len(all_genes), len(sig_genes), float(count)/N, len(lowly_exp_genes), len(undetected_genes)
		# sig_genes_str = ''
		# for g in sig_genes: sig_genes_str += (g + ', ')
		# lowly_exp_genes_str = ''
		# for g in lowly_exp_genes: lowly_exp_genes_str += (g + ', ')
		# undetected_genes_str = ''
		# for g in undetected_genes: undetected_genes_str += (g + ', ')
		# rnaseq_csvwriter.writerow([l] + [sig_genes_str, lowly_exp_genes_str, undetected_genes_str])










# with open(sam_output_file, 'wb') as sam_csvfile:
# 	sam_csvwriter = csv.writer(sam_csvfile)

# with open(rnaseq_output_file, 'wb') as rnaseq_csvfile:
# 	rnaseq_csvwriter = csv.writer(rnaseq_csvfile)

# 	for l in ['NR_034078']:
# 		rnaseq_file = l + '_fold_nofilter_2.csv'
# 		with open(rnaseq_file, 'rb') as gf:
# 			data = [row for row in csv.reader(gf.read().splitlines())]
# 			diff_genes = dict([(row[0], row[1:]) for i,row in enumerate(data[1:])])
# 		baseline_file = l + '_TPM.csv'
# 		with open(baseline_file, 'rb') as gf:
# 			data = [row for row in csv.reader(gf.read().splitlines())]
# 			baseline = dict([(row[0], row[1:]) for i,row in enumerate(data[1:])])
# 		linc_chr, linc_loc = get_TSS(lincs[l])
# 		sig_genes = []
# 		all_genes = []
# 		undetected_genes = []
# 		lowly_exp_genes = []
# 		n_diff = 0
# 		for d in diff_genes.keys():
# 			if l == 'TCONS_00026344':
# 				a = get_activation(diff_genes[d][:7], diff_genes[d][7:])
# 			else:
# 				a = get_activation(diff_genes[d][:9], diff_genes[d][9:])
# 			if a: 
# 				n_diff += 1
# 		for g in genes.keys():
# 			gene_chr, gene_loc = get_TSS(genes[g])

# 			if gene_chr == linc_chr:
# 				if abs(linc_loc - gene_loc) < bp_cutoff:
# 					g = genes_map[g]
# 					if l == 'NR_034078':
# 						print g + '\t' + diff_genes[g]
# 					# if g in sam_hits.keys():
# 					# 	sam_csvwriter.writerow([l, g, genes_map[g], sam_hits[g]])
# 					if g in diff_genes.keys():
# 						# rnaseq_csvwriter.writerow([l, g] + diff_genes[g])
# 						if l == 'TCONS_00026344':
# 							a = get_activation(diff_genes[g][:7], diff_genes[g][7:])
# 						else:
# 							a = get_activation(diff_genes[g][:9], diff_genes[g][9:])
# 						# if a:
# 						if g not in sig_genes: sig_genes.append(g)
# 					else:
# 						if g == 'C6orf195': g = 'LINC01600'
# 						if g == 'FLJ44313': continue
# 						if g == 'LOC101927322': g = 'NM_001289967'
# 						if g == 'CCDC11': g = 'CFAP53'
# 						if g == 'KCNA6': continue
# 						if g == 'GPR133': g = 'ADGRD1'
# 						num_expressed = [0,0,0,0,0]
# 						total_exp = 0
# 						for i, exp in enumerate(baseline[g]):
# 							total_exp += float(exp)
# 							if float(exp) >= 1:
# 								num_expressed[i/3] += 1
# 						if (num_expressed[3] + num_expressed[4]) == 0 and (num_expressed[0] == 0 or num_expressed[1] == 0 or num_expressed[2] == 0):
# 							if g not in undetected_genes and g not in lowly_exp_genes:
# 								if total_exp < 0.001:
# 									undetected_genes.append(g)
# 								else:
# 									lowly_exp_genes.append(g)
# 							# rnaseq_csvwriter.writerow([l, g] + baseline[g] + ['undetected'])
# 					if g not in all_genes: all_genes.append(g)
# 		N = 10000
# 		prob = float(n_diff)/len(genes_map.values())
# 		count = 0
# 		# if l == 'TCONS_00015940': sig_genes = sig_genes + ['MOB3B', 'IFNK', 'EQTN']
# 		for i in range(N):
# 			k = 0
# 			for j in range(len(all_genes)):
# 				if r.random() <= prob: k += 1
# 			if k == len(sig_genes):
# 				count += 1

# 		# if l in rnaseq_lincs:
# 		print l, len(diff_genes.keys()), len(all_genes), len(sig_genes), float(count)/N, len(lowly_exp_genes), len(undetected_genes)
# 		sig_genes_str = ''
# 		for g in sig_genes: sig_genes_str += (g + ', ')
# 		lowly_exp_genes_str = ''
# 		for g in lowly_exp_genes: lowly_exp_genes_str += (g + ', ')
# 		undetected_genes_str = ''
# 		for g in undetected_genes: undetected_genes_str += (g + ', ')
# 		rnaseq_csvwriter.writerow([l] + [sig_genes_str, lowly_exp_genes_str, undetected_genes_str])
