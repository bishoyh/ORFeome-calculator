"""Find ORFs in Genomes and calculate their length."""
from Bio import SeqIO
import matplotlib.pyplot as plt
import sys
import getopt as gtopt



plt.figure(figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xlabel("Orf Size", fontsize=16)
plt.ylabel("Count", fontsize=16)


table = 11
min_pro_len = 60


# Modified From biopython cookbook recipe 19.1.13


def find_orfs_with_trans(seq, trans_table, min_protein_length):
    """Find orfs and translate them.
    :type seq: object
    """

    answer = []
    seq_len = len(seq)

    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    else:
                        start = seq_len - frame - aa_end * 3 - 3
                        end = seq_len - frame - aa_start * 3
                    answer.append((start, end, strand,
                                   trans[aa_start:aa_end]))
                aa_start = aa_end + 1
    answer.sort()
    #print answer
    return answer


def parse(fastafile, output_filename):
    handle = open(fastafile, "rU")  # Change to point to the correct genome file

    L = []
    chars = set('X')
    for record in SeqIO.parse(handle, "fasta"):
        print record.id
        orf_list = find_orfs_with_trans(record.seq, table, min_pro_len)
        for start, end, strand, pro in orf_list:
            if 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' in pro:  # ignore proteins that spanned masked regions
               # print("%s...%s - length %i, strand %i, %i:%i,%s" % (pro[:100], pro[-3:], len(pro) * 3, strand, start, end, record.id))
                print 'Skipped one \n'
                continue
            L.append(len(pro) * 3)
            #print("%s...%s - length %i, strand %i, %i:%i,%s" % (pro[:100], pro[-3:], len(pro) * 3, strand, start, end, record.id))
    handle.close()

    f = open(output_filename, 'w')
    for ele in L:
        f.write(str(ele) + '\n')
    f.close
    plt.hist(L, color="#3F5D7D", bins=100)
    png_filename = output_filename + ".png"
    plt.savefig(png_filename, bbox_inches="tight")


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = gtopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except gtopt.GetoptError:
        print 'ORFer.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'test.py -i <inputfile> -o <outputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    print 'Input file is "', inputfile
    print 'Output file is "', outputfile
    parse(inputfile, outputfile)

if __name__ == "__main__":
    main(sys.argv[1:])
