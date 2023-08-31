from os import close


class DSSPParser:
    file = "ressources/dataset/CATH_info.txt"

    def __init__(self, directory, trainset_size=3000):
        """
        Read the DSSP files and retrieve the amino acid sequences of the chains and their structures
        Create a list of the sequences and structures for the training set and a list for the testing set
        :param directory: str, destination folder where we can find the dssp files
        :param trainset_size: number of proteins in the training set vs the proteins in the testing set
        """

        self.work_dir = directory
        self.trainset_size = trainset_size
        self.trainset = list()
        self.testset = list()
        self.parse()

    def parse_dssp(self, dssp, prot_ID):
        """
        You'll never guess what this method does ;)

        :param dssp: dssp file name
        :param prot_ID: Protein identifier from CATH_info
        :return: (str(amino acid chain), str(structure))
        """

        aa_chain = ""
        structure = ""

        with open(self.work_dir + "/" + dssp + ".dssp", "r") as dssp_file:
            # Ugly checks
            for j, line in enumerate(dssp_file):

                if j < 28:  # Ignore header info
                    continue

                if line[10] != " ":
                    continue

                try:
                    int(line[5:10])
                except ValueError:
                    continue

                if line[13] in "*!BXZ":
                    continue

                if line[11] != prot_ID:
                    continue

                if line[13].islower():
                    aa_chain += "C"
                else:
                    aa_chain += line[13]

                if line[16] in "CTS ":
                    structure += "C"
                elif line[16] in "HGI":
                    structure += "H"
                else:
                    structure += "E"

        return (aa_chain, structure)

    def parse(self):
        """
        Extract the information for each identifiers in the CATH info file. Each identifier contains the id of the
        protein, as well as the chain of the desired sequence. The amino acid X, Z, B will be ignored, and all
        lower caps amino acids are converted to Cystein (C). The structures are converted into C (C, T, S and " "),
        H (H, G, I) or E (E and B).
        :return: None, but create a list of tuples of size 2 (str, str), respectively the amino acid sequence and the
        sequence of the structure
        """

        visited = list()  # Track unique structures
        count = 0

        with open("ressources/dataset/CATH_info.txt", "r") as cath_file:

            for line in cath_file.readlines():
                file_name = line[:4]
                prot_ID = line[4]

                # Done with trainset
                if count == self.trainset_size:
                    self.testset.append(self.parse_dssp(file_name, prot_ID))
                    continue

                # Get trainset data
                self.trainset.append(self.parse_dssp(file_name, prot_ID))
                if file_name not in visited:
                    # Unique structure count
                    visited.append(file_name)
                    count += 1

    def get_trainset(self):
        return self.trainset

    def get_testset(self):
        return self.testset
