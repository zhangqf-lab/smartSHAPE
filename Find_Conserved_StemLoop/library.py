#-*- coding:utf-8 -*-

import Structure, GAP, General, Colors, threading
import getopt, random, os, math, re, sys, time

def load_NCBI_annotation():
    ROOT = "/150T/zhangqf/GenomeAnnotation/NCBI/"
    cattle = GAP.init(ROOT+"cattle.genomeCoor.bed", ROOT+"cattle_transcriptome.fa", showAttr=False)
    chicken = GAP.init(ROOT+"chicken.genomeCoor.bed", ROOT+"chicken_transcriptome.fa", showAttr=False)
    chimpanzee = GAP.init(ROOT+"chimpanzee.genomeCoor.bed", ROOT+"chimpanzee_transcriptome.fa", showAttr=False)
    dog = GAP.init(ROOT+"dog.genomeCoor.bed", ROOT+"dog_transcriptome.fa", showAttr=False)
    frog = GAP.init(ROOT+"frog.genomeCoor.bed", ROOT+"frog_transcriptome.fa", showAttr=False)
    macaque = GAP.init(ROOT+"macaque.genomeCoor.bed", ROOT+"macaque_transcriptome.fa", showAttr=False)
    rat = GAP.init(ROOT+"rat.genomeCoor.bed", ROOT+"rat_transcriptome.fa", showAttr=False)
    zebrafish = GAP.init(ROOT+"zebrafish.genomeCoor.bed", ROOT+"zebrafish_transcriptome.fa", showAttr=False)
    ## More
    # 羊驼
    alpaca = GAP.init(ROOT+"more/alpaca.genomeCoor.bed", ROOT+"more/alpaca_transcriptome.fa", showAttr=False)
    # 狒狒
    baboon = GAP.init(ROOT+"more/baboon.genomeCoor.bed", ROOT+"more/baboon_transcriptome.fa", showAttr=False)
    # 倭黑猩猩
    bonobo = GAP.init(ROOT+"more/bonobo.genomeCoor.bed", ROOT+"more/bonobo_transcriptome.fa", showAttr=False)
    # 几维鸟
    brown_kiwi = GAP.init(ROOT+"more/brown_kiwi.genomeCoor.bed", ROOT+"more/brown_kiwi_transcriptome.fa", showAttr=False)
    # 猫
    cat = GAP.init(ROOT+"more/cat.genomeCoor.bed", ROOT+"more/cat_transcriptome.fa", showAttr=False)
    # 海豚
    dolphin = GAP.init(ROOT+"more/dolphin.genomeCoor.bed", ROOT+"more/dolphin_transcriptome.fa", showAttr=False)
    # 雪貂 
    ferret = GAP.init(ROOT+"more/ferret.genomeCoor.bed", ROOT+"more/ferret_transcriptome.fa", showAttr=False)
    # 乌蛇
    garter_snake = GAP.init(ROOT+"more/garter_snake.genomeCoor.bed", ROOT+"more/garter_snake_transcriptome.fa", showAttr=False)
    # 金雕; 鹫;
    golden_eagle = GAP.init(ROOT+"more/golden_eagle.genomeCoor.bed", ROOT+"more/golden_eagle_transcriptome.fa", showAttr=False)
    # 绿长尾猴
    green_monkey = GAP.init(ROOT+"more/green_monkey.genomeCoor.bed", ROOT+"more/green_monkey_transcriptome.fa", showAttr=False)
    # 狨; 绒
    marmoset = GAP.init(ROOT+"more/marmoset.genomeCoor.bed", ROOT+"more/marmoset_transcriptome.fa", showAttr=False)
    # 裸鼢鼠
    naked_mole_rat = GAP.init(ROOT+"more/naked_mole_rat.genomeCoor.bed", ROOT+"more/naked_mole_rat_transcriptome.fa", showAttr=False)
    # 负鼠
    opossum = GAP.init(ROOT+"more/opossum.genomeCoor.bed", ROOT+"more/opossum_transcriptome.fa", showAttr=False)
    # 鸭嘴兽
    platypus = GAP.init(ROOT+"more/platypus.genomeCoor.bed", ROOT+"more/platypus_transcriptome.fa", showAttr=False)
    # 兔子
    rabbit = GAP.init(ROOT+"more/rabbit.genomeCoor.bed", ROOT+"more/rabbit_transcriptome.fa", showAttr=False)
    # 绵羊
    sheep = GAP.init(ROOT+"more/sheep.genomeCoor.bed", ROOT+"more/sheep_transcriptome.fa", showAttr=False)
    # Ensembl Mouse
    ROOT = "/150T/zhangqf/GenomeAnnotation/Gencode/"
    mouse = GAP.init(ROOT+"mm10.genomeCoor.bed", ROOT+"mm10_transcriptome.fa", showAttr=False)
    human = GAP.init(ROOT+"hg38.genomeCoor.bed", ROOT+"hg38_transcriptome.fa", showAttr=False)
    
    annotation = { 
        "cattle": cattle,
        "chicken": chicken,
        "chimpanzee": chimpanzee,
        "dog": dog,
        "frog": frog,
        "macaque": macaque,
        "human": human,
        "rat": rat,
        "zebrafish": zebrafish,
        "mouse": mouse,
        "alpaca": alpaca,
        "baboon": baboon,
        "bonobo": bonobo,
        "brown_kiwi": brown_kiwi,
        "cat": cat,
        "dolphin": dolphin,
        "ferret": ferret,
        "garter_snake": garter_snake,
        "golden_eagle": golden_eagle,
        "green_monkey": green_monkey,
        "marmoset": marmoset,
        "naked_mole_rat": naked_mole_rat,
        "opossum": opossum,
        "platypus": platypus,
        "rabbit": rabbit,
        "sheep": sheep
    }
    
    return annotation

def load_mouse_gene_description():
    inFile = "/150T/zhangqf/GenomeAnnotation/Gene_Description/mouse.txt"
    geneDescrip = {}
    for line in open(inFile):
        if line.startswith("ENS"):
            data = line.strip().split('\t')
            if len(data) == 3:
                index = data[2].find("[Source:MGI")
                if index != -1:
                    geneDescrip[data[0]] = data[2][:index].strip()
                else:
                    geneDescrip[data[0]] = data[2]
    return geneDescrip

def build_geneSymbol_To_geneID(Parser):
    geneParser = Parser.getGeneParser()
    NameToID = { geneParser[k]['gene_name'].upper(): k for k in geneParser }
    return NameToID

def build_ggt(Homo_Parser):
    for species in Homo_Parser:
        Homo_Parser[species].NameToID = build_geneSymbol_To_geneID(Homo_Parser[species])

def read_MGI_homology(inFile, Parser, Homo_Parser):
    homoTable = {}
    
    geneParser = Parser.getGeneParser()
    NameToID = { geneParser[k]['gene_name']: k for k in geneParser }
    
    hitCount = 0
    loseCount = 0
    lineCount = 0
    noGeneIDCount = { "cattle":0, "chicken":0, "chimpanzee":0, 
                    "dog":0, "frog":0, "macaque":0, "human":0, 
                    "rat":0, "zebrafish":0 }
    
    GeneIDCount = { "cattle":0, "chicken":0, "chimpanzee":0, 
                    "dog":0, "frog":0, "macaque":0, "human":0, 
                    "rat":0, "zebrafish":0 }
    
    for line in open(inFile):
        lineCount += 1
        if lineCount == 1: continue
        
        data = line.strip().split('\t')
        group_id = data[0]
        species_name = data[1].split(',')[0]
        gene_symbol = data[3]
        gene_id = data[4]
        
        if species_name == 'mouse':
            try:
                gene_id = NameToID[gene_symbol]
                gene_types = geneParser[gene_id]['gene_type']
                if 'mRNA' not in gene_types and 'protein_coding' not in gene_types:
                    continue
                hitCount += 1
            except KeyError:
                loseCount += 1
                continue
        else:
            try:
                gene_info = Homo_Parser[species_name].getGeneParser()[gene_id]
                GeneIDCount[species_name] += 1
            except:
                noGeneIDCount[species_name] += 1
                continue
            gene_types = gene_info['gene_type']
            if 'mRNA' not in gene_types and 'protein_coding' not in gene_types:
                continue
        try:
            homoTable[ group_id ].append( (species_name, gene_symbol, gene_id) )
        except KeyError:
            homoTable[ group_id ] = [ (species_name, gene_symbol, gene_id) ]
    
    group_ids = homoTable.keys()
    for group_id in list(group_ids):
        if 'mouse' not in [ it[0] for it in homoTable[group_id] ]:
            del homoTable[group_id]
    
    print(GeneIDCount)
    print(noGeneIDCount)
    print("geneID found: %s;  No geneID found: %s" % (hitCount, loseCount))
    return homoTable

### Get { mouseGeneID => [ (species, geneID), (species, geneID)... ] }
def build_homogene_group(homoTable, Homo_Parser):
    ## 0. Build geneName => geneID dict
    geneName_geneID_dict = build_geneSymbol_To_geneID(Homo_Parser['mouse'])
    
    mouse_dict = {}
    ## 1. Collect mouse genes
    for group_id in homoTable:
        
        ## collect mouse gene names
        mouse_geneNames = []
        for it in homoTable[group_id]:
            species_name, gene_symbol, gene_id = it
            if species_name == 'mouse':
                mouse_geneNames.append(gene_symbol.upper())
        
        ## Skip no mouse gene name items
        if len(mouse_geneNames)==0: continue
        
        ## Generate mouse_gene_id => [ (species_name1, gene_id1), (species_name2, gene_id2), ... ]
        for mouse_geneName in mouse_geneNames:
            try:
                mouse_geneID = geneName_geneID_dict[mouse_geneName]
            except KeyError:
                continue
            for it in homoTable[group_id]:
                species_name, gene_symbol, gene_id = it
                gene_symbol = gene_symbol.upper()
                if species_name == 'mouse':
                    if gene_symbol == mouse_geneName:
                        continue
                    else:
                        try:
                            gene_id = geneName_geneID_dict[gene_symbol]
                        except KeyError:
                            continue
                try:
                    mouse_dict[mouse_geneID].append( (species_name, gene_id) )
                except KeyError:
                    mouse_dict[mouse_geneID] = [ (species_name, gene_id) ]
    
    homogroup = {}
    ## Get all species same name genes
    for gene_id in mouse_dict:
        gene_ensembl = mouse_dict[gene_id]
        ## Get upper gene name
        gene_name = Homo_Parser['mouse'].getGeneParser()[gene_id]['gene_name'].upper()
        for species in Homo_Parser:
            try:
                if species == 'naked_mole_rat' and gene_name == 'Gapdh': print ("here...")
                species_gene_id = Homo_Parser[species].NameToID[gene_name]
                if species == 'naked_mole_rat' and gene_name == 'Gapdh': print (species_gene_id)
            except:
                continue
            gene_types = Homo_Parser[species].getGeneParser()[species_gene_id]['gene_type']
            if 'mRNA' in gene_types or 'protein_coding' in gene_types:
                tuple_in = (species, species_gene_id)
                if tuple_in not in gene_ensembl:
                    gene_ensembl.append( tuple_in )
        homogroup[gene_id] = gene_ensembl
    
    ## Delete those less than 3
    gene_id_list = list(homogroup.keys())
    for gene_id in gene_id_list:
        if len(homogroup[gene_id]) <= 2:
            del homogroup[gene_id]
    
    return homogroup

# ### Get { 'species=tid=geneName':sequence, 'species=tid=geneName':sequence,... }
def build_sequence_db(geneID, Homo_Parser, homogenes, main_species="mouse"):
    homo_species_gid_list = homogenes[geneID] + [ (main_species, geneID) ]
    
    sequence_db = {}
    for species, geneID in homo_species_gid_list:
        trans_list = Homo_Parser[species].getTransByGeneID(geneID)
        for tid in trans_list:
            ft = Homo_Parser[species].getTransFeature(tid)
            if ft['gene_type'] in ('protein_coding', 'mRNA'):
                ID = "%s=%s=%s" % (species, tid, ft['gene_name'])
                try:
                    sequence = Homo_Parser[species].getTransSeq(tid)
                except KeyError:
                    continue
                sequence_db[ID] = sequence
    
    return sequence_db

def trim_matches_for_item(global_matches_item, left_trim, right_trim):
    trimmed_global_matches_item = []
    
    for query_seq, ref_id, species_seq, map_pos in global_matches_item:
        assert len(query_seq) == len(species_seq)
        
        i = 0
        count_base_l = 0
        while count_base_l<left_trim:
            if query_seq[i] != '-':
                count_base_l += 1
            i += 1
        while query_seq[i] == '-':
            i += 1
        
        j = len(query_seq)-1
        count_base_r = 0
        while count_base_r < right_trim:
            if query_seq[j] != '-':
                count_base_r += 1
            j -= 1
        while query_seq[j] == '-':
            j -= 1
        
        trimed_query_seq = query_seq[i:j+1]
        trimed_species_seq = species_seq[i:j+1]
        map_pos += i-species_seq[:i].count('-')
        
        trimmed_global_matches_item.append( [trimed_query_seq, ref_id, trimed_species_seq, map_pos] )
    
    return trimmed_global_matches_item

def trim_matches(global_matches, left_trim, right_trim):
    trimmed_global_matches = {}
    
    for stem_loop in global_matches:
        trimmed_global_matches_item = trim_matches_for_item(global_matches[stem_loop], left_trim, right_trim)
        trimmed_global_matches[stem_loop] = trimmed_global_matches_item
    
    return trimmed_global_matches

def leave_UTR_matches(global_matches, Homo_Parser, UTR='5'):
    filted_matches = {}
    
    for region_str in global_matches:
        filted_matches[region_str] = []
        for mouse_seq, species_transID_geneName, species_seq, map_pos in global_matches[region_str]:
            species, sp_tid, geneName = species_transID_geneName.split('=')
            position = Homo_Parser[species].getRNAPosition(sp_tid, [map_pos, map_pos+len(species_seq)-species_seq.count('-')-1])
            if UTR == '5' and position in ('5UTR', 'span_5UTR_CDS'):
                filted_matches[region_str].append( [mouse_seq, species_transID_geneName, species_seq, map_pos] )
            elif UTR == '3' and position in ('3UTR', 'span_CDS_3UTR'):
                filted_matches[region_str].append( [mouse_seq, species_transID_geneName, species_seq, map_pos] )
            elif UTR == 'CDS' and position in ('span_5UTR_CDS', 'CDS', 'span_CDS_3UTR'):
                filted_matches[region_str].append( [mouse_seq, species_transID_geneName, species_seq, map_pos] )
    
    return filted_matches

def remove_redundancy_for_item(global_matches_item, Homo_Parser):
    filted_global_matches_item = []
    
    ## 1. classify species
    species_dict = {}
    for item in global_matches_item:
        query_seq, species_tid_gName, species_seq, map_pos = item[0], item[1], item[2], int(item[3])
        species, tid, gName = species_tid_gName.split('=')
        
        try:
            species_dict[species].append( [tid, gName, species_seq, map_pos, query_seq] )
        except KeyError:
            species_dict[species] = [ [tid, gName, species_seq, map_pos, query_seq] ]
    
    ## 2. remove redundancy
    for species in species_dict:
        if len(species_dict[species]) == 1:
            tid, gName, species_seq, map_pos, query_seq = species_dict[species][0]
            filted_global_matches_item.append( [query_seq, species+"="+tid+"="+gName, species_seq, map_pos] )
        else:
            seq_dict = {}
            for tid, gName, species_seq, map_pos, query_seq in species_dict[species]:
                try:
                    seq_dict[species_seq].append( [tid, gName, species_seq, map_pos, query_seq] )
                except KeyError:
                    seq_dict[species_seq] = [ [tid, gName, species_seq, map_pos, query_seq] ]
            for species_seq in seq_dict:
                if len(seq_dict[species_seq]) == 1:
                    tid, gName, species_seq, map_pos, query_seq = seq_dict[species_seq][0]
                    filted_global_matches_item.append( [query_seq, species+"="+tid+"="+gName, species_seq, map_pos] )
                else:
                    last_positions = [ ]
                    for tid, gName, species_seq, map_pos, query_seq in seq_dict[species_seq]:
                        chrID, chrPos, Strand = Homo_Parser[species].transCoor2genomeCoor(tid, map_pos)
                        find = False
                        for basicInfo, genomePos in last_positions:
                            if abs(genomePos-chrPos) <= 3:
                                find = True
                                break
                        if not find:
                            item = [[query_seq, species+"="+tid+"="+gName, species_seq, map_pos], chrPos]
                            last_positions.append( item )
                    filted_global_matches_item += [ it[0] for it in last_positions ]
    
    return filted_global_matches_item

def remove_redundancy(global_matches, Homo_Parser):
    filted_global_matches = {}
    
    for stem_loop in global_matches:
        filted_global_matches_item = remove_redundancy_for_item(global_matches[stem_loop], Homo_Parser)
        filted_global_matches[stem_loop] = filted_global_matches_item
    
    return filted_global_matches

def count_covariation(ref_aligned_seq, input_aligned_seq, ref_aligned_dot):
    """
    ref_aligned_seq             -- Reference aligned sequence
    input_aligned_seq           -- An input aligned sequence to annotate
    ref_aligned_dot             -- Reference aligned dotbracket structure
    
    Compare the homologous sequence with reference sequence with known structure
    
    Return a dictionary of:
        - cov: covariation base pairs
        - hf: half flip base pairs
        - bb: breacked base pairs
        - mutl: mutated loop bases
        - insdel: insertion and deletion
        - consbp: conserved base pairs
    """
    
    assert len(ref_aligned_seq) == len(input_aligned_seq) == len(ref_aligned_dot)
    ref_aligned_seq = ref_aligned_seq.replace('U', 'T')
    input_aligned_seq = input_aligned_seq.replace('U', 'T')
    
    anno_align_seq = list(input_aligned_seq)
    
    pair_bps = [ 'AT', 'TA', 'CG', 'GC', 'TG', 'GT' ]
    
    ##counter
    covariation = 0
    half_flip = 0
    breaked_bp = 0
    mut_loop = 0
    insert_del = 0
    conserved_bp = 0
    
    ## Base pairs mutatiom
    bp_list = Structure.dot2ct(ref_aligned_dot)
    for bp in bp_list:
        raw_bp = ref_aligned_seq[ bp[0]-1 ]+ref_aligned_seq[ bp[1]-1 ]
        new_bp = input_aligned_seq[ bp[0]-1 ]+input_aligned_seq[ bp[1]-1 ]
        assert raw_bp in pair_bps
        if new_bp in pair_bps:
            if raw_bp != new_bp:
                if raw_bp[0] == new_bp[0] or raw_bp[1] == new_bp[1]:
                    half_flip += 1
                else:
                    covariation += 1
            else:
                conserved_bp += 1
        else:
            breaked_bp += 1
    
    ## Loop mutation
    for i in range(len(anno_align_seq)):
        if ref_aligned_dot[i] == '.':
            if input_aligned_seq[i] != ref_aligned_seq[i]:
                if input_aligned_seq[i] != '-' and ref_aligned_seq[i] != '-':
                    mut_loop += 1
    
    ## insertion/Deletion
    for b1,b2 in zip(ref_aligned_seq, input_aligned_seq):
        if (b1 == '-' and b2 != '-') or (b1 != '-' and b2 == '-'):
            insert_del += 1
    
    return { 'cov':covariation, 'hf':half_flip, 'bb':breaked_bp, 'mutl':mut_loop, 'insdel':insert_del, 'consbp':conserved_bp}

def calc_conserve_score(ref_aligned_seq, input_aligned_seq, ref_aligned_dot):
    """
    ref_aligned_seq             -- Reference aligned sequence
    input_aligned_seq           -- An input aligned sequence to annotate
    ref_aligned_dot             -- Reference aligned dotbracket structure
    
    Calculate three scores:
        break_ratio             -- Ratio of breaked base paires
        conserved_bp_ratio      -- Ratio of conserved base pairs
        covariation_ratio       -- Ratio of covariate bases in paired bases
    
    Return (break_ratio, conserved_bp_ratio, covariation_ratio)
    """
    
    total_bases = len(ref_aligned_seq)
    
    cov_dict = count_covariation(ref_aligned_seq, input_aligned_seq, ref_aligned_dot)
    paired_num = (cov_dict['cov']+cov_dict['hf']+cov_dict['bb']+cov_dict['consbp'])*2
    break_ratio = (2.0*cov_dict['bb']+cov_dict['insdel'])/paired_num
    bp_cons_ratio = 2.0*(cov_dict['cov']+cov_dict['hf']+cov_dict['consbp'])/paired_num
    cov_ratio = (2.0*cov_dict['cov']+1.0*cov_dict['hf'])/(2.0*cov_dict['consbp']+2.0*cov_dict['cov']+2.0*cov_dict['hf']+0.01)
    
    return round(break_ratio,3), round(bp_cons_ratio,3), round(cov_ratio,3)

def calc_pred_cons_structure_score(ref_aligned_seq, input_aligned_seq, ref_aligned_dot, max_loop_bias=2):
    """
    ref_aligned_seq             -- Reference aligned sequence
    input_aligned_seq           -- An input aligned sequence to annotate
    ref_aligned_dot             -- Reference aligned dotbracket structure
    max_loop_bias               -- Maximun center loop bias
    
    Calculate three scores:
        loop_is_in_frame            -- True or False
        compact_ratio               -- compact_input/compact_ref
        structure_simi_ratio        -- Ratio of co-paired bases/all paired bases
        diff_loop_size              -- Input loop size - reference loop size
        pred_aligned_structure      -- Dotbracket structure
    
    Return (loop_is_in_frame, compact_ratio, structure_simi_ratio, diff_loop_size, pred_aligned_structure)
    """
    
    assert len(ref_aligned_seq) == len(input_aligned_seq) == len(ref_aligned_dot)
    clean_input_seq = input_aligned_seq.replace("-", "")
    dot = Structure.predict_structure(clean_input_seq)
    if '(' not in dot:
        return False, 0, 0, 0, "."*len(input_aligned_seq)
    
    aligned_dot = Structure.dot_to_alignDot(dot, input_aligned_seq)
    
    slim_aligned_dot = ""
    slim_ref_aligned_dot = ""
    for s1,s2 in zip(aligned_dot, ref_aligned_dot):
        if s1!='-' or s2!='-':
            slim_aligned_dot += s1
            slim_ref_aligned_dot += s2
    
    ## Reference loop
    ref_loop_start = slim_ref_aligned_dot.rindex("(")+2
    ref_loop_end = slim_ref_aligned_dot.index(")")
    ## Input loop
    inp_loop_start = slim_aligned_dot.rindex("(")+2
    inp_loop_end = slim_aligned_dot.index(")")
    
    #######
    ### For loop size calculation
    #######
    ref_loop_size = slim_ref_aligned_dot[ref_loop_start-1:ref_loop_end].count('.')
    input_loop_size = slim_aligned_dot[inp_loop_start-1:inp_loop_end].count('.')
    diff_loop_size = input_loop_size - ref_loop_size
    
    ## Loop is in frame ??
    if abs(ref_loop_start-inp_loop_start)<=max_loop_bias and abs(ref_loop_end-inp_loop_end)<=max_loop_bias:
        loop_in_frame = True
    else:
        loop_in_frame = False
    
    ## Compact ratio
    tmp_clean_aligned_dot = slim_ref_aligned_dot.replace("-", "")
    ref_compact = 2.0*tmp_clean_aligned_dot.count("(")/len(tmp_clean_aligned_dot)
    input_compact = 2.0*dot.count("(")/len(dot)
    compact_ratio = round(input_compact/ref_compact,3)
    
    ## Structure similarity
    total = 0
    same_paired = 0
    for s1,s2 in zip(slim_aligned_dot, slim_ref_aligned_dot):
        if s1 == '-' and s2  == '-':
            continue
        if s1=='(' or s2 == ')':
            if s1==s2:
                same_paired += 1
            total += 1
    
    structure_simi_ratio = round(1.0*same_paired/total,3)
    
    #print slim_aligned_dot
    #print slim_ref_aligned_dot
    
    return loop_in_frame, compact_ratio, structure_simi_ratio, diff_loop_size, aligned_dot

def filter_bad_consist_structure_for_item(global_matches_item, ref_raw_dot, Homo_Parser):
    """
    global_matches_item             -- { '1864-1882' => [ ['ACTGGAAAGTAACTC-CAGT', 'rat=XM_006231248.2=Cd274', 'ACTGCTAAGCAAGTTACCCA', 1146]... ] }
    ref_raw_dot                     -- Dot bracket structure of raw reference sequence
    Homo_Parser                     -- { 'mouse':mouse_gaper, 'human':human_gaper, 'frog':frog_gaper,... }
    
    """
    
    new_list = []
    for ref_aligned_seq,input_id,input_aligned_seq,input_pos in global_matches_item:
        
        #print ref_raw_dot, ref_aligned_seq
        ref_aligned_dot = Structure.dot_to_alignDot( ref_raw_dot, ref_aligned_seq )
        bb_ratio, bp_cons_ratio, cov_ratio = calc_conserve_score(ref_aligned_seq, input_aligned_seq, ref_aligned_dot)
        if bb_ratio > 0.2: continue
        
        loop_is_in_frame, compact_ratio, structure_simi_ratio, diff_loop_size, pred_slim_aligned_dot = calc_pred_cons_structure_score(ref_aligned_seq, input_aligned_seq, ref_aligned_dot, max_loop_bias=2)
        if not loop_is_in_frame: continue
        if compact_ratio < 0.7: continue
        if structure_simi_ratio < 0.7: continue
        if diff_loop_size >= 2: continue
        if pred_slim_aligned_dot.count('(') < 5: continue
        
        pure_structure = pred_slim_aligned_dot.replace('-', '')
        new_list.append( [ref_aligned_seq,input_id,input_aligned_seq,input_pos,pure_structure] )
    
    return new_list

def filter_bad_consist_structure(global_matches, structure_dict, Homo_Parser, min_align_num=2):
    filted_good_matches = {}
    
    for stem_loop in global_matches:
        #print stem_loop
        good_list = filter_bad_consist_structure_for_item(global_matches[stem_loop], structure_dict[stem_loop], Homo_Parser)
        if len(good_list) >= min_align_num:
            filted_good_matches[stem_loop] = good_list
    
    return filted_good_matches

def unify_search_alignment(global_matches):
    unified_alignments = {}
    
    for query_id in global_matches:
        #print query_id
        
        alignments = global_matches[query_id]
        align_num = len(alignments)
        correct_ref = [""] * align_num
        positions = [0] * align_num
        
        if align_num <= 1: continue
        
        new_query = ""
        pure_query = alignments[0][0].replace("-","")
        
        while len(new_query.replace("-",""))!=len(pure_query):
            bases = []
            
            for i in range(len(positions)):
                bases.append( alignments[i][0][positions[i]] )
            if 1<=bases.count('-')<len(bases):
                for i,base in enumerate(bases):
                    if base == '-':
                        correct_ref[i] += alignments[i][2][positions[i]]
                        positions[i] += 1
                    else:
                        correct_ref[i] += '-'
                new_query += '-'
            else:
                uniq_base = list(set(bases))
                assert len(uniq_base) == 1
                new_query += uniq_base[0]
                for i in range(len(bases)):
                    correct_ref[i] += alignments[i][2][positions[i]]
                    positions[i] += 1
        
        unified_alignments[query_id] = []
        for i in range(len(alignments)):
            unified_alignments[query_id].append( [new_query, alignments[i][1], correct_ref[i], alignments[i][3]]+alignments[i][4:] )
    
    return unified_alignments

def output_fine_alignment(transID, unified_alignments, structure_dict, shape_dict, Homo_Parser, geneDescription, main_species='mouse', writePredict=False, OUT=sys.stdout):
    for region_str in sorted(unified_alignments.keys(), key=lambda x: int(x.split('-')[0])):
        ls, re = region_str.split('-')
        ls, re = int(ls), int(re)
        ref_dot = structure_dict[region_str]
        ref_shape = shape_dict[region_str]
        
        ref_aligned_seq = unified_alignments[region_str][0][0]
        ref_aligned_dot = Structure.dot_to_alignDot(ref_dot, ref_aligned_seq)
        ref_aligned_shape = Structure.shape_to_alignSHAPE(ref_shape, ref_aligned_seq)
        
        geneID = Homo_Parser[main_species].getTransFeature(transID)['gene_id'].split('.')[0]
        score = unified_alignments[region_str][0][-1]
        if geneID in geneDescription:
            print(">%s\t%s\t%s\t%s" % (transID, region_str, geneDescription[geneID], score), file=OUT)
        else:
            print(">%s\t%s\t%s" % (transID, region_str, score), file=OUT)
        seq_struc_pairs = []
        for data in unified_alignments[region_str]:
            assert data[0] == ref_aligned_seq
            input_id, input_aligned_seq, input_pos = data[1], data[2], data[3]
            
            input_annoned_seq = Structure.annotate_covariation(ref_aligned_seq, input_aligned_seq, ref_aligned_dot, anno_loop=True)
            input_sp,input_transID,input_geneName = input_id.split('=')
            if main_species in input_id: input_id = Colors.f(input_id, 'red')
            
            position_label_str = Homo_Parser[input_sp].labelRNAPosition(input_transID,[input_pos, input_pos])
            format_string = "%s\t%s\t%s\n" % (input_annoned_seq, position_label_str, input_id)
            OUT.writelines(format_string)
            
            if len(data) >= 5:
                input_dot = data[4].replace("-", "")
                input_aligned_dot = Structure.dot_to_alignDot(input_dot, input_aligned_seq)
                seq_struc_pairs.append( (input_id, input_aligned_seq, input_aligned_dot) )
        
        print(ref_aligned_dot, file=OUT)
        print(Colors.color_SHAPE(ref_aligned_shape), file=OUT)
        print("", file=OUT)
        
        if writePredict:
            for input_id, input_aligned_seq, input_aligned_dot in seq_struc_pairs:
                print(input_aligned_seq+'\t'+input_id, file=OUT)
                print(input_aligned_dot, file=OUT)
            print("", file=OUT)
        
        if len(unified_alignments[region_str]) > 4:
            smoothed_conserv, position_str = get_conserve_profile( unified_alignments[region_str], Homo_Parser )
            print("## Sequence conservation profile", file=OUT)
            print(Colors.color_SHAPE(smoothed_conserv, cutoff=[0.7, 0.8, 0.9]), file=OUT)
            print(position_str, file=OUT)
            print("", file=OUT)

def calc_aligned_conserv_score(aligned_seq_list, binNum=50, binSize=False):
    import numpy as np
    
    seqNum = len(aligned_seq_list)
    seqLen = len(aligned_seq_list[0])
    for seq in aligned_seq_list[1:]:
        assert len(seq) == seqLen
    
    ## produce nucleotide-based conservation
    conserv_list = []
    for i in range(seqLen):
        base_list = [ seq[i] for seq in aligned_seq_list ]
        base_list.sort()
        A = base_list.count('A')
        T = base_list.count('T')
        C = base_list.count('C')
        G = base_list.count('G')
        Gap = base_list.count('-')
        A /= float(seqNum)
        T /= float(seqNum)
        C /= float(seqNum)
        G /= float(seqNum)
        Gap /= float(seqNum)
        conserv = max([A, T, C, G])
        assert A+T+C+G != 0.0
        if Gap > 0.5:
            conserv_list.append("NULL")
        else:
            conserv_list.append( conserv )
    
    if binSize:
        wSize = binSize
        binNum = int(round(1.0*seqLen/wSize))
    else:
        assert seqLen > binNum
        wSize = 1.0*seqLen/binNum
        binNum = int(binNum)
    
    i = 0
    smoothed_conserv = []
    while len(smoothed_conserv)<binNum:
        start = int(round(i))
        end = int(round(i+wSize))
        sub_conv = conserv_list[start:end]
        if sub_conv.count('NULL') / len(sub_conv) >= 0.5:
            smoothed_conserv.append( 'NULL' )
        else:
            sub_conv = [ it for it in sub_conv if it != 'NULL' ]
            smoothed_conserv.append( round(np.mean(sub_conv),3) )
        i += wSize
    
    return smoothed_conserv

def get_conserve_profile( unified_alignments_item, Homo_Parser ):
    full_seq_list = []
    subseq_list = []
    i = 0
    while i < len(unified_alignments_item):
        ref_seq,sp_tid,sp_seq,sp_pos,sp_ss = unified_alignments_item[i][:5]
        species,tid,gene_name = sp_tid.split('=')
        trans_seq = Homo_Parser[species].getTransSeq(tid)
        ft = Homo_Parser[species].getTransFeature(tid)
        
        sp_seq = sp_seq.replace("-", "")
        start = trans_seq.find(sp_seq)
        end = start + len(sp_seq)
        
        full_seq_list.append( trans_seq[max(start-50,0):end+50] )
        subseq_list.append(sp_seq)
        
        i += 1
    
    aligned_seq_list = Structure.multi_alignment(full_seq_list)
    smoothed_conserv = calc_aligned_conserv_score(aligned_seq_list, binNum=50)
    bin_len = len(aligned_seq_list[0])/50.0
    start, end = Structure.align_find(aligned_seq_list[0], subseq_list[0])
    start_bin = int(start/bin_len)
    end_bin = int(end/bin_len)
    position_str = " "*start_bin + "-"*(end_bin-start_bin+1) + " "*(49-end_bin)
    return smoothed_conserv, position_str

def isYRY(subseq):
    if re.match(".*[TC][AG][TC]", subseq):
        return True
    return False

def process_UTR(tid, shape_list, seq_db, Homo_Parser, main_species='mouse', UTR='5', min_structure_score=0.65,
    min_identity=0.6, ext=3, ERR=sys.stderr):
    
    assert UTR in ('5', '3', 'CDS')
    
    sequence = Homo_Parser[main_species].getTransSeq(tid)
    ft = Homo_Parser[main_species].getTransFeature(tid)
    geneID =  ft['gene_id']
    cds_start, cds_end = ft['cds_start'], ft['cds_end']
    tLen = ft['trans_len']
    
    if UTR == '5':
        no_null = cds_start - shape_list[:cds_start].count('NULL')
        if cds_start<100:
            print("Error: %s has a too short 5'UTR" % (geneID, ), file=ERR)
            return [], {}, {}, {}
        if no_null<40:
            print("Error: %s has too few shape scores in 5'UTR" % (geneID, ), file=ERR)
            return [], {}, {}, {}
    elif UTR == '3':
        no_null = tLen - cds_end -  shape_list[cds_end:].count('NULL')
        if tLen-cds_end<100:
            print("Error: %s has a too short 3'UTR" % (geneID, ), file=ERR)
            return [], {}, {}, {}
        if no_null<40:
            print("Error: %s has too few shape scores in 3'UTR" % (geneID, ), file=ERR)
            return [], {}, {}, {}
    elif UTR == 'CDS':
        no_null = cds_end - cds_start -  shape_list[cds_start:cds_end].count('NULL')
        if cds_end-cds_start<100:
            print("Error: %s has a too short CDS" % (geneID, ), file=ERR)
            return [], {}, {}, {}
        if no_null<40:
            print("Error: %s has too few shape scores in CDS" % (geneID, ), file=ERR)
            return [], {}, {}, {}
    
    if UTR == '5':
        stem_loops = Structure.sliding_score_stemloop(sequence, shape_list, start=0, end=cds_start, wSize=100, wStep=50, skip_noshape=True)
    elif UTR == '3':
        stem_loops = Structure.sliding_score_stemloop(sequence, shape_list, start=cds_end, wSize=100, wStep=50, skip_noshape=True)
    else:
        stem_loops = Structure.sliding_score_stemloop(sequence, shape_list, start=cds_start, end=cds_end, wSize=100, wStep=50, skip_noshape=True)
    
    stem_loops = [ it for it in stem_loops if it[2] >= min_structure_score ]
    
    ## filter YRY
    i = 0
    while i<len(stem_loops):
        r,dot,score = stem_loops[i]
        subseq = sequence[r[1]-1:r[2]]
        if isYRY(subseq):
            i += 1
        else:
            del stem_loops[i]
    
    #print( "len(stem_loops)", len(stem_loops), stem_loops )
    if len(stem_loops) == 0:
        print("Error: %s -- no stem-loop found" % (geneID, ), file=ERR)
        return [], {}, {}, {}
       
    UTR_seq_db = {}
    for seq_id in seq_db:
        species, sp_id, sp_gene_name = seq_id.split('=')
        sp_ft = Homo_Parser[species].getTransFeature(sp_id)
        sp_cds_start = sp_ft['cds_start']
        sp_cds_end = sp_ft['cds_end']
        sp_len = sp_ft['trans_len']
        if UTR=='5' and sp_cds_start > 30:
            UTR_seq_db[seq_id] = seq_db[seq_id]
        elif UTR=='3' and sp_len-sp_cds_end > 30:
            UTR_seq_db[seq_id] = seq_db[seq_id]
        else:
            UTR_seq_db[seq_id] = seq_db[seq_id]
    if len(UTR_seq_db) <= 3:
        print("Error: %s has too few homologous 5'UTR sequence %s" % (geneID, len(UTR_seq_db)), file=ERR)
        return [], {}, {}, {}
    
    ####### Organize stem loop
    sl_motif_dict = {}      # stem-loop motif
    sl_structure_dict = {}  # stem-loop structure
    sl_shape_dict = {}      # stem-loop shape
    
    for region,ss,score in stem_loops:
        ls,le,rs,re = region
        if ls-1-ext>=0 and re+ext<=tLen:
            key = "%s-%s" % (ls, re)
            sl_motif_dict[ key ] = sequence[ls-1-ext:re+ext]
            sl_structure_dict[ key ] = ss
            sl_shape_dict[ key ] = shape_list[ls-1:re]
    
    ####### Search homologous
    global_matches = Structure.global_search(sl_motif_dict, UTR_seq_db, evalue=20, min_identity=min_identity, thread_nums=3)
    trimmed_global_matches = trim_matches(global_matches, left_trim=ext, right_trim=ext)
    filt_5UTR_global_matches = leave_UTR_matches(trimmed_global_matches, Homo_Parser, UTR=UTR)
    filt_global_matches =  remove_redundancy(filt_5UTR_global_matches, Homo_Parser)
    fine_global_matches = filter_bad_consist_structure(filt_global_matches, sl_structure_dict, Homo_Parser, min_align_num=3)
    unified_alignments = unify_search_alignment(fine_global_matches)
    
    #if len(unified_alignments) < 4:
    #    print("Error: %s has too few homologous structures" % (geneID, ), file=ERR)
    #    return [], {}, {}, {}
    #print( "len(stem_loops)", len(stem_loops), stem_loops )
    for region,ss,score in stem_loops:
        ls,le,rs,re = region
        if ls-1-ext>=0 and re+ext<=tLen:
            key = "%s-%s" % (ls, re)
            if score>0.9 and key not in unified_alignments:
                gene_name = Homo_Parser[main_species].getTransFeature(tid)['gene_name']
                ID = f'{main_species}={tid}={gene_name}'
                s = sequence[ls-1:re]
                unified_alignments[key] = [ [s, ID, s, ls, ss, score] ]
                sl_motif_dict[ key ] = sequence[ls-1-ext:re+ext]
                sl_shape_dict[ key ] = shape_list[ls-1:re]
                sl_structure_dict[ key ] = ss
            elif key in unified_alignments:
                unified_alignments[key][0].append(score)
    
    #print(unified_alignments)
    #print( "unified_alignments", unified_alignments.keys() )
    return unified_alignments, sl_motif_dict, sl_structure_dict, sl_shape_dict

def find_stable_conserved_stemloop_for_trans(transID, Homo_Parser, geneDescription, shape_dict, homogenes, 
    main_species='mouse', find_region=['5', '3'], OUT=sys.stdout, 
    ERR=sys.stderr, min_structure_score=0.65, min_identity=0.6, 
    global_lock=None, ext=3):
    """
    transID                         -- Transcript ID
    Homo_Parser                     -- Multiple homologous GAPers
    shape_dict                      -- SHAPE dictionary
    mouseGeneIDToSpeciesInfo        -- object returned by build_mouse_gene_group
    min_identity                    -- Minimun identify for global_search
    global_lock                     -- For multiple-thread process
    """
    ft = Homo_Parser[main_species].getTransFeature(transID)
    geneID = ft['gene_id']
    if geneID not in homogenes:
        print("Error: %s not in homogenes" % (geneID, ), file=ERR)
        return []
    
    if ft['gene_type'] not in ('mRNA', 'protein_coding'):
        print("Error: %s is not a mRNA" % (geneID, ), file=ERR)
        return []
    
    cds_start = ft['cds_start']
    cds_end = ft['cds_end']
    tLen = ft['trans_len']
    
    seq_db = build_sequence_db(geneID, Homo_Parser, homogenes, main_species=main_species)
    if len(seq_db) <= 3:
        print("Error: %s has too few homologous sequence %s" % (geneID, len(seq_db)), file=ERR)
        return []
    
    shape_list = shape_dict[transID]
    sequence = Homo_Parser[main_species].getTransSeq(transID)
    if len(shape_list) != tLen or len(sequence) != tLen:
        print("Error: %s have different length -- annoLen:%s; shapeLen:%s; seqLen:%s" % (geneID, tLen, len(shape_list), len(sequence)), file=ERR)
        return []
    
    stem_loop_found = {}
    for region in find_region:
        unified_alignments, sl_motif_dict, sl_structure_dict, sl_shape_dict = process_UTR(transID, shape_list, seq_db, Homo_Parser, main_species=main_species, UTR=region, min_identity=min_identity, ext=ext, min_structure_score=min_structure_score, ERR=ERR)
        stem_loop_found[region] = [unified_alignments, sl_motif_dict, sl_structure_dict, sl_shape_dict]
    
    #print "predict 5"
    #unified_alignments_5, sl_motif_dict_5, sl_structure_dict_5, sl_shape_dict_5 = process_UTR(transID, shape_list, seq_db, Homo_Parser, main_species=main_species, UTR='5', min_identity=min_identity, ext=ext, min_structure_score=min_structure_score, ERR=ERR)
    #print "predict 3"
    #unified_alignments_3, sl_motif_dict_3, sl_structure_dict_3, sl_shape_dict_3 = process_UTR(transID, shape_list, seq_db, Homo_Parser, main_species=main_species, UTR='3', min_identity=min_identity, ext=ext, min_structure_score=min_structure_score, ERR=ERR)
    
    if global_lock:
        while global_lock.locked():
            continue
        global_lock.acquire()
    
    #print( "len(stem_loop_found)", len(stem_loop_found) )
    for region in stem_loop_found:
        unified_alignments, sl_motif_dict, sl_structure_dict, sl_shape_dict = stem_loop_found[region]
        if len(unified_alignments) >= 1:
            output_fine_alignment(transID, unified_alignments, sl_structure_dict, sl_shape_dict, Homo_Parser, geneDescription, main_species=main_species, writePredict=True, OUT=OUT)
    
    if global_lock:
        global_lock.release()
    
    return unified_alignments


#### Load Annoation
Homo_Parser = load_NCBI_annotation()
geneDescription = load_mouse_gene_description()

#### Build geneSymbol_To_geneID for each species
build_ggt(Homo_Parser)

#### Load homologous genes
homoFile = "/150T/zhangqf/GenomeAnnotation/homology/MGI/HOM_AllOrganism.rpt"
MGI_homoTable = read_MGI_homology(homoFile, Homo_Parser['mouse'], Homo_Parser)
homogenes = build_homogene_group(MGI_homoTable, Homo_Parser)


"""
for transID in set(mouse_smartSHAPE):
    stemloop = find_stable_conserved_stemloop_for_trans(transID, Homo_Parser, geneDescription, mouse_smartSHAPE, homogenes, main_species='mouse',
                                                    OUT=sys.stdout, ERR=sys.stderr, min_structure_score=0.65,
                                                    min_identity=0.6, global_lock=None, ext=3)
"""

