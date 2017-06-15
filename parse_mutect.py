import sys


def parse_mutect(mutect_output):
    """ Parse mutect output 
    """

    mutect_hash = {}
    with open(mutect_output,"r") as IN:
        for line in IN:
            line=line.strip('\n')
            if line[0] == "#":
                continue
            contents=line.split("\t")
            if contents[0] == "contig": ## get the columns
                mutect_columns = contents
            
            else:
                key = contents[0]+"_"+contents[1]
                mutect_hash[key] = line

    return [mutect_columns,mutect_hash]
                
def parse_smcounter(smcounter_output):
    """ Parse smcounter output
    """
    smcounter_hash = {}
    with open(smcounter_output,"r") as IN:
        for line in IN:
            line=line.strip('\n')
            contents = line.split("\t")
            if contents[0] == "CHROM":
                smcounter_columns = contents
            else:
                key = contents[0]+"_"+contents[1]
                smcounter_hash[key] = line
                
    return [smcounter_columns,smcounter_hash]

def parse_falsepositives(falsepositives_vcf):
    """ Parse the false positives vcf file
    """

    false_positives_hash = {}
    with open(falsepositives_vcf,"r") as IN:
        for line in IN:
            line=line.strip('\n')
            if line[0] == "#":
                continue
            else:
                contents = line.split('\t')
                key = contents[0]+"_"+contents[1]
                false_positives_hash[key] = line

    return false_positives_hash

def join_mutect_smcounter(falsepositives_vcf,mutect_output,smcounter_output):
    """ Join smcounter and mutect results
    """
    
    false_positives_hash = parse_falsepositives(falsepositives_vcf)
    mutect_columns,mutect_hash = parse_mutect(mutect_output)
    smcounter_columns,smcounter_hash = parse_smcounter(smcounter_output)
    combined_columns = mutect_columns+smcounter_columns
    #print "\t".join(combined_columns)
    ref = 0
    alt = 0
    count = 0
    #print smcounter_columns
    print mutect_columns
    for loci in false_positives_hash:
        if loci in smcounter_hash:
            smcounter_contents = smcounter_hash[loci].split('\t')
            mutect_contents =  mutect_hash[loci].split('\t')
            alt_count = float(mutect_contents[mutect_columns.index("t_alt_count")])
            ref_count = float(mutect_contents[mutect_columns.index("t_ref_count")])
            mutect_af = float(mutect_contents[mutect_columns.index("tumor_f")])
            smcounter_af = float(smcounter_contents[smcounter_columns.index("sVMF")])
            dp = alt_count+ref_count
            #print mutect_af,smcounter_af,dp,alt_count
            if mutect_af < 0.01:
                count+=1
            
            if smcounter_contents[2] != mutect_contents[3]:
                ref+=1
                
            if smcounter_contents[3] != mutect_contents[4]: ## Alts dont match 
                alt+=1
            print smcounter_contents[3],mutect_contents[4],mutect_contents[2],mutect_contents[mutect_columns.index("t_ref_max_mapq")],mutect_contents[mutect_columns.index("t_alt_max_mapq")],mutect_contents[-1],mutect_contents[-2]
            combined_contents = mutect_hash[loci]+'\t'+smcounter_hash[loci]
        else:
            print "Zz"
                
            #print combined_contents
    #print count
    

if len(sys.argv) != 4:
    print "RUN as : python parse_mutect.py <false_positives.vcf> <mutect_output> <smcounter_output>\n"
    sys.exit(-1)

join_mutect_smcounter(sys.argv[1],sys.argv[2],sys.argv[3])
            
