import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QTextEdit, QComboBox, QPushButton, QFileDialog
from PyQt5.QtGui import QColor, QTextCursor
from PyQt5.QtCore import Qt
from collections import Counter

class DNAProcessor(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()

        self.DNA_Codons = {
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "TGT": "C", "TGC": "C",
            "GAT": "D", "GAC": "D",
            "GAA": "E", "GAG": "E",
            "TTT": "F", "TTC": "F",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
            "CAT": "H", "CAC": "H",
            "ATA": "I", "ATT": "I", "ATC": "I",
            "AAA": "K", "AAG": "K",
            "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "ATG": "M",
            "AAT": "N", "AAC": "N",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAA": "Q", "CAG": "Q",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "TGG": "W",
            "TAT": "Y", "TAC": "Y",
            "TAA": "_", "TAG": "_", "TGA": "_"
        }

        self.RNA_Codons = {
            "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "UGU": "C", "UGC": "C",
            "GAU": "D", "GAC": "D",
            "GAA": "E", "GAG": "E",
            "UUU": "F", "UUC": "F",
            "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
            "CAU": "H", "CAC": "H",
            "AUA": "I", "AUU": "I", "AUC": "I",
            "AAA": "K", "AAG": "K",
            "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
            "AUG": "M",
            "AAU": "N", "AAC": "N",
            "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAA": "Q", "CAG": "Q",
            "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
            "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
            "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
            "UGG": "W",
            "UAU": "Y", "UAC": "Y",
            "UAA": "_", "UAG": "_", "UGA": "_"
        }

        self.AminoAcidNames = {
            "A": "Alanine", "C": "Cysteine", "D": "Aspartic Acid", "E": "Glutamic Acid",
            "F": "Phenylalanine", "G": "Glycine", "H": "Histidine", "I": "Isoleucine",
            "K": "Lysine", "L": "Leucine", "M": "Methionine", "N": "Asparagine",
            "P": "Proline", "Q": "Glutamine", "R": "Arginine", "S": "Serine",
            "T": "Threonine", "V": "Valine", "W": "Tryptophan", "Y": "Tyrosine",
            "_": "Stop"
        }

        self.Codon_table = {**self.DNA_Codons, **self.RNA_Codons}
    

    def initUI(self):
        self.setWindowTitle("Sequence Analysis")
        self.setGeometry(300, 300, 400, 300)

        layout = QVBoxLayout()
        self.setLayout(layout)

        label = QLabel("Nucleotide Sequence Analysis Tool")
        layout.addWidget(label)

        self.sequenceEdit = QLineEdit()
        self.sequenceEdit.setPlaceholderText("Enter your DNA or RNA sequence")
        layout.addWidget(self.sequenceEdit)

        self.choiceCombo = QComboBox()
        self.choiceCombo.addItem("Sequence Length")
        self.choiceCombo.addItem("Nucleotide Frequency")
        self.choiceCombo.addItem("Convert to RNA")
        self.choiceCombo.addItem("GC Content")
        self.choiceCombo.addItem("Reading Frames")
        self.choiceCombo.addItem("Protein")
        self.choiceCombo.addItem("Triplet Analysis")
        layout.addWidget(self.choiceCombo)

        self.resultText = QTextEdit()
        self.resultText.setReadOnly(True)
        layout.addWidget(self.resultText)

        button = QPushButton("Run")
        button.clicked.connect(self.run)
        layout.addWidget(button)

        exit_button = QPushButton("Exit", self)
        exit_button.clicked.connect(self.close)
        layout.addWidget(exit_button)

        self.show()

    def close(self):
        self.close()

    def run(self):
        sequence = self.sequenceEdit.text().upper()
        choice = self.choiceCombo.currentText()

        if not sequence:
            self.resultText.setText("Please enter a sequence")
            return

        sequence_type = self.validate_sequence(sequence)
        if not sequence_type:
            self.resultText.setText("Invalid sequence. Please ensure all nucleotides are A, C, G, T, or U and in uppercase.")
            return

        is_rna = sequence_type == 'RNA'

        if choice == "Sequence Length":
            self.resultText.setText(f"Sequence Length: {len(sequence)}")

        elif choice == "Nucleotide Frequency":
            frequency = {nucleotide: sequence.count(nucleotide) for nucleotide in set(sequence)}
            self.resultText.setText(f"Nucleotide Frequency: {frequency}")

        elif choice == "Convert to RNA":
            if not is_rna:
                rna_sequence = sequence.replace('T', 'U')
                self.resultText.setText(f"RNA Sequence: {rna_sequence}")
            else:
                self.resultText.setText("The sequence is already RNA.")

        elif choice == "GC Content":
            gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
            self.resultText.setText(f"GC Content: {gc_content:.2f}%")

        elif choice == "Reading Frames":
            frames = self.reading_frames(sequence)
            self.resultText.setText("Reading frames for sequence:")
            for i, frame in enumerate(frames, 1):
                self.resultText.append(f"Frame {i}: {(''.join(frame)) if frame else 'None'}")

        elif choice == "Protein":
            proteins = self.get_proteins(sequence)
            self.resultText.setText("Proteins:")
            for protein in proteins:
                self.resultText.append(protein)

        elif choice == "Triplet Analysis":
            triplets = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
            self.resultText.setText(f"Triplets: {triplets}")

        elif choice == "Amino Acid Sequence":
            amino_acid_sequence = " ".join(self.translate_seq(sequence, 0))
            self.resultText.setText(f"Amino Acid sequence from DNA: {amino_acid_sequence}")
            amino_acid = input("Enter the CODON to get the frequency: ")
            codon_frequency = self.codon_usage(sequence, amino_acid)
            self.resultText.append(f"Codon frequency for {amino_acid}: {codon_frequency}")

    def validate_sequence(self, sequence):
        if 'T' in sequence:
            print("Sequence is DNA")
            return 'DNA'
        elif 'U' in sequence:
            print("Sequence is RNA")
            return 'RNA'
        print("Invalid sequence. Please ensure all nucleotides are A, C, G, T, or U and in uppercase.")
        return None

    def DNART(self, seq):
        seq = seq.upper()  # Convert to uppercase
        if 'U' in seq:
            # RNA complement
            complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
        else:
            # DNA complement
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join([complement[nuc] for nuc in seq])[::-1]

    def translate_seq(self, seq, init_pos=0):
        # Make sure the sequence length is a multiple of 3
        remainder = len(seq) % 3
        if remainder!= 0:
            seq = seq + 'N' * (3 - remainder)  # 'N' denotes any nucleotide, filling to make the length a multiple of 3

        # Translate the sequence up to the last complete codon
        return [self.Codon_table.get(seq[pos:pos + 3], 'X') for pos in range(init_pos, len(seq), 3)]

    def codon_usage(self, seq, aminoacid):
        tmpList = []
        for i in range(0, len(seq) - 2, 3):
            if self.Codon_table.get(seq[i:i + 3]) == aminoacid:
                tmpList.append(seq[i:i + 3])

        freqDict = dict(Counter(tmpList))
        totalweight = sum(freqDict.values())
        result = {}
        for s in freqDict:
            result[s] = (self.AminoAcidNames[self.Codon_table[s]], round(freqDict[s] / totalweight, 2))   
        return result

    def reading_frames(self, seq):
        frames = []
        frames.append(self.translate_seq(seq, 0))
        frames.append(self.translate_seq(seq, 1))
        frames.append(self.translate_seq(seq, 2))
        frames.append(self.translate_seq(self.DNART(seq), 0))
        frames.append(self.translate_seq(self.DNART(seq), 1))
        frames.append(self.translate_seq(self.DNART(seq), 2))
        return frames

    def get_proteins(self, sequence):
        sequence = str(sequence)  # Convert sequence to string if it's not already
        frames = self.reading_frames(sequence.upper())
        proteins = []
        for frame in frames:
            proteins.extend(self.proteins_rf(frame))
        return proteins

    def proteins_rf(self, aa_seq):
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def all_proteins_from_orfs(self, seq, startReadPos=0, endReadPos=0, ordered=False):
        if endReadPos > startReadPos:
            rfs = self.reading_frames(seq[startReadPos: endReadPos])
        else:
            rfs = self.reading_frames(seq)

        res = []
        for rf in rfs:
            prots = self.proteins_rf(rf)
            for p in prots:
                res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)
        return res

if __name__ == "__main__":
    app = QApplication(sys.argv)
    dna_processor = DNAProcessor()
    sys.exit(app.exec_())