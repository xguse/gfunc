from spartan.utils.files import ParseFastA,fold_seq


def make_fasta_from_gene_names_with_Tx_or_Prot_dbs(gene_names,out_path,fasta_path_list=[]):
    """
    """
    # create one large nested fasta dict keyed as fastas[gene_name][transcript_prot_name]=seq
    
    fasta_db = dev.dict_tree()
    
    for fasta_path in fasta_path_list:
        p = ParseFastA(fasta_path)
        d = p.to_dict()
        for name,seq in d.iteritems():
            fasta_db[name[:-3]][name] = seq
    
    # collect my seqs
    my_seqs = dev.dict_tree()
    
    for name in gene_names:
        children = fasta_db[name]
        for child in children.keys():
            my_seqs[child] = fasta_db[name][child]
            
    # write seqs to file
    write_fasta_from_seq_dict(my_seqs,out_path)
    
    
def write_fasta_from_seq_dict(seq_dict,out_path):
    """
    """
    out_file = open(out_path,'w')
    
    for name,seq in seq_dict.iteritems():
        folded = fold_seq(seq)
        folded = '\n'.join(folded)
        rec = '>%s\n%s\n' % (name,folded)
        out_file.write(rec)
        
    out_file.close()


def write_pos_neg_seq_fastas(pos_names,fasta_paths,pos_path,neg_path,neg_names='left_overs',write_negs=True):
    """
    function to output forground/background fasta files for motif discovery purposes
    """
    
    
    # build mega-seq-dict
    mega_dict = {}
    
    for fasta_path in fasta_paths:
        p = ParseFastA(fasta_path)
        d = p.to_dict()
        mega_dict.update(d)
    
    # if neg_names = left_overs retrieve neg_name list
    if neg_names == 'left_overs':
        neg_names = list(set(mega_dict.keys()) - set(pos_names))
        
    # build pos_dict
    pos_dict = {}
    for name in pos_names:
        pos_dict[name] = mega_dict[name]

    # build neg_dict
    if write_negs:
        neg_dict = {}
        for name in neg_names:
            neg_dict[name] = mega_dict[name]
        
    # write out_files
    write_fasta_from_seq_dict(pos_dict,pos_path)
    if write_negs:
        write_fasta_from_seq_dict(neg_dict,neg_path)    
    