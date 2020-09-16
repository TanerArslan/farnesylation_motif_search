from Bio import SeqIO
'''
Specific motifs were sarched based on input motif.
C: Cysteine
X: Any amino acid but C
'''

class MotifSearch(object):

    def __init__(self, filename):
        self.filename = filename
        self.fastaSequences = SeqIO.parse(open(self.filename), 'fasta')

    # Motif for CCXXX
    def ccxxx(self):

        self.ccxxxfile = open('Motif_CCxxx_amino_N_terminal.txt', 'w')

        for self.fasta in self.fastaSequences:
            name, seq, des = self.fasta.id, str(self.fasta.seq), str(self.fasta.description)
            if len(seq) > 5 and seq[-5] == "C" and seq[-4] == "C" and "C" not in seq[-3:]:
                self.ccxxxfile.write(
                    name.split('|')[1] + "\t" + des.partition("GN=")[2].split(" ")[0] + "\t" + seq[-5:] + "\n")
            else:
                pass

        self.ccxxxfile.close()

        return()

    # Motif for CXCXX
    def cxcxx(self):

        self.cxcxxfile = open('Motif_CxCxx_amino_N_terminal.txt', 'w')

        for self.fasta in self.fastaSequences:
            name, seq, des = self.fasta.id, str(self.fasta.seq), str(self.fasta.description)
            if len(seq) > 5 and seq[-5] == "C" and seq[-3] == "C" and "C" not in seq[-4] and "C" not in seq[-2:]:
                self.cxcxxfile.write(
                    name.split('|')[1] + "\t" + des.partition("GN=")[2].split(" ")[0] + "\t" + seq[-5:] + "\n")
            else:
                pass

        self.cxcxxfile.close()

        return()

    # Motif for CXXCX
    def cxxcx(self):

        self.cxxcxfile = open('Motif_CxxCx_amino_N_terminal.txt', 'w')

        for self.fasta in self.fastaSequences:
            name, seq, des = self.fasta.id, str(self.fasta.seq), str(self.fasta.description)
            if len(seq) > 5 and seq[-5] == "C" and seq[-2] == "C" and "C" not in seq[-4:-2] and "C" not in seq[-1:]:
                self.cxxcxfile.write(
                    name.split('|')[1] + "\t" + des.partition("GN=")[2].split(" ")[0] + "\t" + seq[-5:] + "\n")
            else:
                pass

        self.cxxcxfile.close()

        return()

    # Motif for CXXXC
    def cxxxc(self):

        self.cxxxcfile = open('Motif_CxxxC_amino_N_terminal.txt', 'w')

        for self.fasta in self.fastaSequences:
            name, seq, des = self.fasta.id, str(self.fasta.seq), str(self.fasta.description)
            if len(seq) > 5 and seq[-5] == "C" and seq[-1] == "C" and "C" not in seq[-4:-1]:
                self.cxxxcfile.write(
                    name.split('|')[1] + "\t" + des.partition("GN=")[2].split(" ")[0] + "\t" + seq[-5:] + "\n")
            else:
                pass

        self.cxxxcfile.close()

        return()

    # Motif for XCCXX
    def xccxx(self):

        self.xccxxfile = open('Motif_xCCxx_amino_N_terminal.txt', 'w')

        for self.fasta in self.fastaSequences:
            name, seq, des = self.fasta.id, str(self.fasta.seq), str(self.fasta.description)
            if len(seq) > 5 and "C" not in seq[-5] and seq[-4:-2] == "CC" and "C" not in seq[-2:]:
                self.xccxxfile.write(
                    name.split('|')[1] + "\t" + des.partition("GN=")[2].split(" ")[0] + "\t" + seq[-5:] + "\n")
            else:
                pass

        self.xccxxfile.close()

        return()

    # Motif for XCXCX
    def xcxcx(self):

        self.xcxcxfile = open('Motif_xCxCx_amino_N_terminal.txt', 'w')

        for self.fasta in self.fastaSequences:
            name, seq, des = self.fasta.id, str(self.fasta.seq), str(self.fasta.description)
            if len(seq) > 5 and "C" not in seq[-5] and seq[-4] == "C" and "C" not in seq[-3] and seq[-2] == "C" and "C" not in seq[-1]:
                self.xcxcxfile.write(
                    name.split('|')[1] + "\t" + des.partition("GN=")[2].split(" ")[0] + "\t" + seq[-5:] + "\n")
            else:
                pass

        self.xcxcxfile.close()

        return()

    # Motif for XCXXC
    def xcxxc(self):

        self.xcxxcfile = open('Motif_xCxxC_amino_N_terminal.txt', 'w')

        for self.fasta in self.fastaSequences:
            name, seq, des = self.fasta.id, str(self.fasta.seq), str(self.fasta.description)
            if len(seq) > 5 and "C" not in seq[-5] and seq[-4] == "C" and "C" not in seq[-3:-1] and seq[-1] == "C":
                self.xcxxcfile.write(
                    name.split('|')[1] + "\t" + des.partition("GN=")[2].split(" ")[0] + "\t" + seq[-5:] + "\n")
            else:
                pass

        self.xcxxcfile.close()

        return()

    # Motif for CCXXX
    def xxccx(self):

        self.xxccxfile = open('Motif_xxCCx_amino_N_terminal.txt', 'w')

        for self.fasta in self.fastaSequences:
            name, seq, des = self.fasta.id, str(self.fasta.seq), str(self.fasta.description)
            if len(seq) > 5 and "C" not in seq[-5:-3] and seq[-3:-1] == "CC" and "C" not in seq[-1]:
                self.xxccxfile.write(
                    name.split('|')[1] + "\t" + des.partition("GN=")[2].split(" ")[0] + "\t" + seq[-5:] + "\n")
            else:
                pass

        self.xxccxfile.close()

        return()

    # Motif for XXCXC
    def xxcxc(self):

        self.xxcxcfile = open('Motif_xxCxC_amino_N_terminal.txt', 'w')

        for self.fasta in self.fastaSequences:
            name, seq, des = self.fasta.id, str(self.fasta.seq), str(self.fasta.description)
            if len(seq) > 5 and "C" not in seq[-5:-3] and seq[-3] == "C" and "C" not in seq[-2] and seq[-1] == "C":
                self.xxcxcfile.write(
                    name.split('|')[1] + "\t" + des.partition("GN=")[2].split(" ")[0] + "\t" + seq[-5:] + "\n")
            else:
                pass

        self.xxcxcfile.close()

        return()

    # Motif for XXXCC
    def xxxcc(self):

        self.xxxccfile = open('Motif_xxxCC_amino_N_terminal.txt', 'w')

        for self.fasta in self.fastaSequences:
            name, seq, des = self.fasta.id, str(self.fasta.seq), str(self.fasta.description)
            if len(seq) > 5 and "C" not in seq[-5:-2] and seq[-2:] == "CC":
                self.xxxccfile.write(
                    name.split('|')[1] + "\t" + des.partition("GN=")[2].split(" ")[0] + "\t" + seq[-5:] + "\n")
            else:
                pass

        self.xxxccfile.close()

        return()


m = MotifSearch('uniprotproteomeAUP000005640.fasta')
#f1 = m.ccxxx()
#f2 = m.cxcxx()
#f3 = m.cxxcx()
#f4 = m.cxxxc()
#f5 = m.xccxx()
#f6 = m.xcxcx()
#f7 = m.xcxxc()
#f8 = m.xxccx()
#f9 = m.xxcxc()
#f10 = m.xxxcc()
