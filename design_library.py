#!/usr/bin/env python
"""
simple of set of python scripts to get target guides
"""
import argparse, json, tempfile, os,itertools, subprocess
import conf_linc
import twobitreader
import sqlite3
import numpy
import time
import Bio.SeqIO
import math
from operator import itemgetter, attrgetter
import csv
# import matplotlib.pyplot as plt

# define parameters
GC_cutoff = 25 #pick guides with GC content % above this value
settings = conf_linc.get_settings()
initial_spacing = 20 #no overlap is tolerated initially
left_overlap = ("tttcttggctttatatatcttGTGGAAAGGACGAAACACC").upper()
right_overlap = ("gttttagagctaggccAACATGAGGATCACC").upper()

#Helper functions
def revcomp(sequence):
    '''
    returns the reverse complement of sequence
    '''
    basecomplement = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'} 
    letters = list(sequence) 
    letters.reverse() 
    dna = ''
    for base in letters:
        dna += basecomplement[base] 
    return dna 

def indexList(s, item, i=0):
    '''
    make a list if indexes of item in s
    '''
    i_list = []
    while True:
        try:
            i = s.index(item, i)
            i_list.append(i)
            i += 1
        except:
            break
    return i_list

def Target_stretch(guide):
    
    st1 = ('AAAA')
    st2 = ('TTTT')
    st3 = ('GGGG')
    st4 = ('CCCC')

    if not ((st1 in guide) or (st2 in guide) or (st3 in guide) or (st4 in guide)):
        return 'true'


def GC_content(guide):
    """
    takes guide sequence as imput, returns true if GC content above
    threshold defined above
    """
    N = guide.count("G")
    N += guide.count("C")
    percent = float(N)/20*100
    if percent > GC_cutoff:
        return 'true'
    else: 
        return 'false'

def get_b_guides(s, ind):
    """
    takes a sequence s and and a list of indices in sequence (indicating start of 'CC')
    returns a list of 20bp bottom guide sequences (with Target_stretch and GC content criteria)
    """
    guides = []
    ind_ver = []
    for i in ind:
        guide = (revcomp(s[i+3:i+23]))
        if len(guide) == 20 and 'N' not in guide:
            if  GC_content(guide) == 'true' and Target_stretch(guide) == 'true':
                 ind_ver.append(i)
                 guides.append(guide) 
    return ind_ver, guides 

def get_t_guides(s, ind):
    """
    takes a sequence s and and a list of indices in sequence (indicating start of 'GG')
    returns a list of 20bp top strand guide sequences (with Target_stretch and GC content criteria)
    """
    guides = []
    ind_ver = []
    for i in ind:
        if i > 20:
            guide = (s[i-21:i-1])
            if len(guide) == 20 and 'N' not in guide:
                if  GC_content(guide) == 'true' and Target_stretch(guide) == 'true':
                    ind_ver.append(i)
                    guides.append(guide) 
    return ind_ver, guides 

def get_guides(seq):
    '''
    finds all top and bottom guides in seq and returns them with indices
    '''
    i_list_CC = indexList(seq, 'CC') #find all the indexes of CC's
    i_list_GG = indexList(seq, 'GG') #find all the indexes of GG's
    (b_guide_i, b_guides) = get_b_guides(seq, i_list_CC)  
    (t_guide_i, t_guides) = get_t_guides(seq, i_list_GG)
    return b_guide_i, b_guides, t_guide_i, t_guides

def get_distance(strand, b_guide_i, t_guide_i, window, bin_offset):
    '''
    returns lists of bottom guide and top guide distances from TSS
    '''
    if strand == "+":
        if bin_offset > 0:
            if bin_offset != 800:
                b_guide_dist = [(bin_offset + 23 + window)-(x+23) for x in b_guide_i]
                t_guide_dist = [(bin_offset + 23 + window)-(x-1) for x in t_guide_i]
            else:
                b_guide_dist = [(bin_offset + window)-(x+23) for x in b_guide_i]
                t_guide_dist = [(bin_offset + window)-(x-1) for x in t_guide_i]
        else:
            b_guide_dist = [(settings.UPSTREAM_REGION_SIZE + window)-(x+23) for x in b_guide_i]
            t_guide_dist = [(settings.UPSTREAM_REGION_SIZE + window)-(x-1) for x in t_guide_i]

    elif strand == "-":
        b_guide_dist = [bin_offset+x+3 for x in b_guide_i]
        t_guide_dist = [bin_offset+x-21 for x in t_guide_i]
    return b_guide_dist, t_guide_dist

def add_location(strand, all_guides, gene):
    '''
    adds the location of each guide in the genome
    '''
    for guide in all_guides:
        if strand == "+":
            guide.extend([gene["chrom"],long(gene["txStart"]) - guide[0]])
        elif strand == "-":
            guide.extend([gene["chrom"],long(gene["txEnd"]) + guide[0]])
        

def get_sorted_guides(seq, gene, window, spacing, bin_offset):
    '''
    returns a list with the remaining guides in seq
    format of lists:
    [distance, sequence, strand (b/t)]
    '''
    all_b_guides = []
    all_t_guides = []
    spaced_guides = []
    guide_count = 0 
    distance = 0
    strand = gene["strand"]

    (b_guide_i, b_guides, t_guide_i, t_guides) = get_guides(seq)
    (b_guide_dist, t_guide_dist) = get_distance(strand, b_guide_i, t_guide_i, window, bin_offset)
    for i in range (0,len(b_guide_dist)): 
        all_b_guides.append([b_guide_dist[i],b_guides[i],"b"])  #makes a nested list of bottom spacers with: [distance, sequence, strand (b/t)]
    for i in range (0,len(t_guide_dist)):
        all_t_guides.append([t_guide_dist[i],t_guides[i],"t"])

    all_guides = all_b_guides + all_t_guides #makes one nested list of bottom and top guides
    add_location(strand, all_guides, gene)
    all_guides_sorted = sorted(all_guides, key=itemgetter(0))  #sorts by distance

    #remove overlapping guides
    for guide in all_guides_sorted:
        guide_count += 1
        if guide_count == 1:
            spaced_guides.append(guide)
            distance = guide[0]
        elif guide[0] > (distance + spacing):
            spaced_guides.append(guide)
            distance = guide[0]
    return spaced_guides

def get_ot_guides(guide_list):
    '''
    connects to db to fetch ot scores for guides in guide_list
    returns a sorted list of guides (according to OT score)
    '''
    guides_scored = []
    sequences = []

    for guide in guide_list:
        # guides_scored.append([guide[1],0.943]) #just for testing
        sequences.append(guide[1])

    dbpath = os.path.abspath(settings.DATABASE_FILE)
    dbaddress = "//{0}".format(dbpath)
    conn = sqlite3.connect(dbaddress)
    c = conn.cursor()
    query = "select * FROM sgrnas WHERE seq IN ({0}) ORDER BY score".format(','.join(['?']*len(sequences)))
    try:
        c.execute(query, sequences)
        guides_scored = c.fetchall() 
       
    except:
        print "did not find guides"

    return guides_scored

def get_bin(guide, bins):
    for loc in bins:
        if guide[0] <= loc:
            return bins.index(loc)

def bin_guides(guides, bins):
    binned_guides = [[] for i in range(len(bins))]
    for guide in guides:
        binned_guides[get_bin(guide,bins)].append(guide)
    binned_count = [len(bin) for bin in binned_guides]
    return binned_guides, binned_count

def bin_region(region, bins, strand):
    binned_region = []
    if strand == "+":
        for i in range(len(bins)):
            if i == 0:
                binned_region.append(region[:bins[i]])
            elif i == len(bins)-1:
                binned_region.append(region[bins[i-1]-20:])
            else:
                binned_region.append(region[bins[i-1]-20:bins[i]])
        binned_region.reverse()
    elif strand == "-":
        for i in range(len(bins)):
            if i == 0:
                binned_region.append(region[:bins[i]+20])
            elif i == len(bins)-1:
                binned_region.append(region[bins[i]:])
            else:
                binned_region.append(region[bins[i-1]:bins[i]+20])
    return binned_region

def add_distances(temp_ot_guides, guide_locations):
    chrom = temp_ot_guides[0][5]
    for guide in temp_ot_guides:
        if chrom in guide_locations.keys():
            guide_locations[chrom].append(guide[-1])
        else:
            guide_locations[chrom] = []

def get_location_bins(bins, gene):
    if gene["strand"] == "+":
        TSS = long(gene["txStart"])-1
        return [TSS - bin for bin in bins]
    elif gene["strand"] == "-":
        TSS = long(gene["txEnd"])
        return [bin + TSS for bin in bins]

def get_adjustment(location_bins, locations,strand):
    adjustment = [0 for bin in location_bins]
    if strand == "+":
        for loc in locations:
            if loc <= location_bins[0] + 200:
                for bin in location_bins:
                    if loc >= bin:
                        index = location_bins.index(bin)
                        adjustment[index] += 1
                        break
    elif strand == "-":
        for loc in locations:
            if loc >= location_bins[0] - 200:
                for bin in location_bins:
                    if loc <= bin:
                        index = location_bins.index(bin)
                        adjustment[index] += 1
                        break
    return adjustment


def adjust_bins(ideal_bins, bins, guide_locations, gene):
    location_bins = get_location_bins(bins,gene)
    if gene["chrom"] not in guide_locations:
        return ideal_bins
    else:
        adjustment = get_adjustment(location_bins, guide_locations[gene["chrom"]], gene["strand"])
        adjusted_bins = []
        for i in range(len(ideal_bins)):
            diff = ideal_bins[i] - adjustment[i]
            if diff >= 0:
                adjusted_bins.append(diff)
            else:
                adjusted_bins.append(0)
        return adjusted_bins

def filter_distances(ot_guides, guide_locations):
    filtered_guides = []
    for guide in ot_guides:
        chrom = guide[5]
        loc = guide[6]
        if chrom not in guide_locations.keys():
            return ot_guides
        elif loc not in guide_locations[chrom]:
            filtered_guides.append(guide)
    return filtered_guides

def list_sgrnas(genes_file, genes_format, genome_2bit):
    '''
    Returns a list of (ontarget) sgrna sequences using a genome file and list of
    transcription start sites form a .csv file.
    '''
    tbf = twobitreader.TwoBitFile(genome_2bit)
    window_ext_count = 0
    isoform_count = 0
    final_ot_guides = []
    all_guides = [] # for calculating ot
    guide_locations = {}
    b_count = 0
    t_count = 0
    overlap_count = 0
    initial_guide_count = []

    #bins for distributing sgRNA positions
    #numbers reflect the larger edge of the respective bin
    #ideal situation: 6 in first bin, 1 in each of the other bins
    bins = [200, 300, 400, 600, 800]
    num_ten = 0
    num_notten = 0

    with open(genes_file, 'rb') as gf:
        f = [row for row in csv.reader(gf.read().splitlines())]
        for i,l in enumerate(f): # i is index, l is entry
            if i == 0: 
                columns = l
                continue
            if settings.MAX_SGRNAS > 0:
                if i > settings.MAX_SGRNAS: break

            #fetch the current gene
            gene = dict([(columns[i],e) for i,e in enumerate(l)])
            isoform_count += 1
            window = 0 #if there are not enough guides in the initial window, this tracks window extension (bp)
            ot_guides = [[] for bin in bins]


            #use strand to compute TSS upstream region bounds
            if gene["strand"] == "+":
                region_bounds = [long(gene["txStart"])-1-settings.UPSTREAM_REGION_SIZE, 
                                     long(gene["txStart"])-1]
            elif gene["strand"] == "-":
                region_bounds = [long(gene["txEnd"]), 
                                     long(gene["txEnd"]) + settings.UPSTREAM_REGION_SIZE]

            #extract upstream region from the genome
            region = tbf[gene["chrom"]][region_bounds[0]:region_bounds[1]]
            region = region.upper()
            #print region
            if "N" in region: 
                print "found N in target region of", gene["name"]
            #     continue

            ideal_bins_std = [6,1,1,1,1]
            ideal_bins = adjust_bins(ideal_bins_std, bins, guide_locations, gene)
            #get >3 guides
            (remaining_guides) = get_sorted_guides(region, gene, window, initial_spacing, 0)

            #return a list of list of guides binned according to bins
            binned_guides, binned_count = bin_guides(remaining_guides, bins)
            binned_region = bin_region(region, bins, gene["strand"])
            current_guides = binned_guides[0]
            current_region = binned_region[0]
            #deal with cases where there are not enough guides in first 200bp (<3)
            #keep increasing the overlap allowed up to 10bp max
            overlap = 0 
            trial = 0

            while len(current_guides) < 3 and trial < 6:
                trial += 1 
                overlap += 2
                if gene["strand"] == "+": bin_offset = bins[0]
                else: bin_offset = 0
                current_guides = get_sorted_guides(current_region, gene, window, initial_spacing - overlap, bin_offset)
            binned_guides[0] = current_guides
            binned_count[0] = len(current_guides)
            #deal with cases where there are not enough guides in first 800bp (<10)
            bin = 0

            while numpy.sum(binned_count) < numpy.sum(ideal_bins):
                if bin < len(bins):
                    current_guides = binned_guides[bin]
                    current_region = binned_region[bin]
                    trial = 0
                    overlap = 0
                    while numpy.sum(binned_count) < numpy.sum(ideal_bins) and trial < 6:
                        trial += 1
                        overlap += 2
                        if bin > 0:
                            if gene["strand"] == "+": bin_offset = bins[bin]
                            else: bin_offset = bins[bin-1]
                        else:
                            if gene["strand"] == "+": bin_offset = bins[0]
                            else: bin_offset = 0
                        current_guides = get_sorted_guides(current_region, gene, window, initial_spacing - overlap, bin_offset)
                        binned_guides[bin] = current_guides
                        binned_count[bin] = len(current_guides)
                    bin += 1
                else: #if neccessary, extend the search region
                    window += 50
                    if gene["strand"] == "+":
                        region_bounds = [long(gene["txStart"])-1-(settings.UPSTREAM_REGION_SIZE + window), 
                                            long(gene["txStart"])-1]
                    elif gene["strand"] == "-":
                        region_bounds = [long(gene["txEnd"]), 
                                            long(gene["txEnd"]) + (settings.UPSTREAM_REGION_SIZE + window)]
                    region = tbf[gene["chrom"]][region_bounds[0]:region_bounds[1]]
                    region = region.upper()
                    current_guides = get_sorted_guides(current_region, gene, window, initial_spacing - overlap, bins[-1])
                    binned_guides[-1] = current_guides
                    binned_count[-1] = len(current_guides)
                if window > 200:
                    print "not enough guides found", gene["name"]
                    print binned_count
                    break

            if 0 < trial < 6:
                overlap_count += 1
            elif window > 0:
                window_ext_count += 1
            
            temp_ot_guides = []

            for bin in range(len(bins)):

                remaining_guides = binned_guides[bin]
                #fetch ot_guides        
                ot_guides_sql = get_ot_guides(remaining_guides) #format: [(u'TTGGGAATCCTGAGTCCAAG', 0.9430855725845865) etc]
                for guide in ot_guides_sql:
                    for r_guide in remaining_guides:
                        if r_guide[1] == str(guide[0]):
                            r_guide.insert(0,guide[1])
                            r_guide.insert(0,gene["name"])
                            ot_guides[bin].append(r_guide)
                ot_guides[bin] = filter_distances(ot_guides[bin], guide_locations)
                ot_guides[bin] = sorted(ot_guides[bin], key=itemgetter(1), reverse=True)
                all_guides.extend(ot_guides[bin])
                #check if enough guides are present in the bin of ot_guides to fulfill the ideal number of guides
                if len(ot_guides[bin]) < ideal_bins[bin]:
                    temp_ot_guides.extend(ot_guides[bin])
                    #if not enough guides in current bin, take more guides in the next bin
                    if bin < len(bins) - 1:
                        ideal_bins_std[bin+1] += (ideal_bins[bin]-len(ot_guides[bin]))
                        ideal_bins_std[bin] = len(ot_guides[bin])
                        ideal_bins = adjust_bins(ideal_bins_std, bins, guide_locations, gene)
                    #if no guides in last bin, take guide from previous bin if available
                    else: 
                        guides_needed = ideal_bins[bin] - len(ot_guides[bin])
                        while guides_needed > 0 and bin >= 0:
                            ideal_bin = ideal_bins[bin]
                            if len(ot_guides[bin]) > ideal_bin:
                                diff = len(ot_guides[bin]) - ideal_bin
                                if diff > guides_needed:
                                    diff = guides_needed
                                temp_ot_guides.extend(ot_guides[bin][ideal_bin:ideal_bin+diff])
                                ideal_bins[bin] += 1
                                guides_needed -= diff
                            bin -= 1                          
                else:
                    temp_ot_guides.extend(ot_guides[bin][:ideal_bins[bin]])
            if len(temp_ot_guides) > 0: add_distances(temp_ot_guides, guide_locations)
            final_ot_guides.extend(temp_ot_guides)
            if len(temp_ot_guides) != 10: num_notten += 1
            else: num_ten += 1

        print num_notten
        print num_ten
        print overlap_count
        print window_ext_count

        print len(final_ot_guides)
        all_guides_tmp = []
        for guide in all_guides:
            all_guides_tmp.append(">dummy")
            all_guides_tmp.append(guide[3])
        print len(all_guides_tmp)
        return final_ot_guides, all_guides_tmp
            

            # if isoform_count % 1000 == 0:
            #     print isoform_count, "isoforms processed"
            #     print window_ext_count, "isoforms with window extension"
            #     print overlap_count, "isoforms with higher overlap"

            # #initial_guide_count.append(len(remaining_guides))
                
    """
    #write resulting ot guides into csv file
    result_file = "web_guides_with_scores.csv"
    with open(result_file,'w') as csvfile:
        mywriter = csv.writer(csvfile)
        for record in final_ot_guides:
            mywriter.writerow(record)
            

    plt.hist(initial_guide_count, bins=8)
    plt.title("Histogram")
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.show()

    """

def writecsv(data, filename):
    with open(filename, 'wb') as csvfile:
        csvwriter = csv.writer(csvfile)
        for row in data:
            csvwriter.writerow(row)      

def __main__():
    print "computing a list of sgrnas"
    print "genes from: {0}, format {1}"\
        .format(settings.GENES_FILE, settings.GENES_FORMAT)
    print "genome: {0}".format(settings.GENOME_2BIT_FILE)
    final_ot_guides, all_guides = list_sgrnas(settings.GENES_FILE, settings.GENES_FORMAT, settings.GENOME_2BIT_FILE)
    final_ot_guides = sorted(final_ot_guides,key=itemgetter(-2,-1,0))
    # all_guides = sorted(all_guides,key=itemgetter(-2,-1,0))
    writecsv(final_ot_guides,'final_linc_guides_v2.csv')
    # writecsv(all_guides, 'all_linc_guides.csv')
    print left_overlap
    print right_overlap

if __name__ == "__main__":
    __main__()

