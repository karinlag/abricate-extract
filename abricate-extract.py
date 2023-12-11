import os, sys, argparse
from Bio import SeqIO, Seq

class specificGene:
    "This class defines a specific gene, mainly taken from abricate output"
    def __init__(self, line, stripFromGenomeName):
        splitline = line.split()
        self.genome = splitline[0]
        self.contig = splitline[1]
        self.start = int(splitline[2])
        self.end = int(splitline[3])
        self.strand = splitline[4]
        self.gene = splitline[5]
        self.coverage = float(splitline[9])
        self.identity = float(splitline[10])
        self.database = splitline[11]
        self.accno = splitline[12]

        # I want to put the name of the genome in the fasta heading, in a nice way
        self.genomeName = os.path.basename(self.genome).rstrip(stripFromGenomeName)
        
        # Here I create the basic seqrecord that will be used for output. 
        self.geneDescription = "__".join([self.contig, self.gene, self.database, self.accno, str(self.start), \
                                    str(self.end), str(self.coverage), str(self.identity)])
        self.id = "__".join([self.gene, self.genomeName])
        self.seqRecord = SeqIO.SeqRecord(id = self.id, description = self.geneDescription, seq= "")

def getGenes(abricateFile, geneToExtract, stripFromGenomeName):
    """
    This method takes in an abricate file, a text string specifying which gene to extract,
    and a file ending to strip from the genome name to get nice genome names.
    """
    openedFile = open(abricateFile, "r")
    lines = openedFile.readlines()
    if "FILE" in lines[0]:
        lines.pop(0) #skip headings
    genes = set()
    for line in lines:
        newGene = specificGene(line, stripFromGenomeName)
        if newGene.gene == geneToExtract:
            genes.add(newGene)  
    return genes

def extractGenes(allMatchingGenes):
    """
    This method takes all genes from the input, grabs the genome file and extracts
    the gene that we are looking for. The sequence is reverse complemented if needed
    and shoved into the Seq part of the seqRecord of the gene.
    """
    extractedGenes = []
    for gene in allMatchingGenes:
        openedFile = open(gene.genome, "r")
        for seq_record in SeqIO.parse(gene.genome, "fasta"):
            if seq_record.name == gene.contig:
                seqToExtract = seq_record.seq[gene.start-1:gene.end]
                if gene.strand == "-":
                    seqToExtract = seqToExtract.reverse_complement()
                gene.seqRecord.seq=seqToExtract
                extractedGenes.append(gene.seqRecord)
    return extractedGenes


def main(abricateFile, geneToExtract, stripFromGenomeName,outputFile):
    allMatchingGenes = getGenes(abricateFile, geneToExtract, stripFromGenomeName)
    foundGeneList = extractGenes(allMatchingGenes)
    SeqIO.write(foundGeneList, outputFile, "fasta")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Extract sequences from abricate output. Fasta files must be reachable as mentioned in abricate output.')
    parser.add_argument('abricateFile', help="abricate file name")   
    parser.add_argument('geneToExtract', help="gene to extract from abricate output, column = GENE")   
    parser.add_argument('stripFromGenomeName', help="text to be removed from genome basename")   
    parser.add_argument('outputFile', help="output file name")   
    args = parser.parse_args()

    main(args.abricateFile, args.geneToExtract, args.stripFromGenomeName, args.outputFile)