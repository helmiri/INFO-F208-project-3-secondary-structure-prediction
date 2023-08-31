from math import log, log10, sqrt
import numpy as np


class GOR:
    def __init__(self, Parser):
        """
        GOR III implementation, divided in a training part and a prediction. The training creates the
        counters/frequency used afterwards to predict which structure an amino acid is part of according to the
        propensity of this amino acid R to be part of a structure S (I(S,R)) and its neighbourhood if it is part
        of this structure (I(S,R,Rj) with -8<=j<=8).
        :param Parser: DSSPParser object
        """
        self.counters = dict()   # {"counter": value,...}
        self.parser = Parser
        for seq_struct in Parser.get_trainset():
            self.update_counters(seq_struct[0], seq_struct[1])

    def individual_inf(self, amino_acid, structure):
        structures = ["C", "H", "E"]
        structures.remove(structure)

        key = "{},{}"

        f_nSR = self.counters[key.format(
            structures[0], amino_acid)] + self.counters[key.format(structures[1], amino_acid)]

        f_nS = self.counters[structures[0]] + self.counters[structures[1]]
        f_S = self.counters[structure]
        f_SR = self.counters[key.format(structure, amino_acid)]

        return log10(f_SR/f_nSR) + log10(f_nS/f_S)

    def neighborhood_inf(self, amino_acid, structure, neighbor, position):
        structures = ["C", "H", "E"]
        structures.remove(structure)

        key1 = "{},{}"
        key2 = "{},{},{},{}"
        f_SRR = self.counters[key2.format(
            structure, neighbor, amino_acid, position)]

        f_SR = self.counters[key1.format(structure, amino_acid)]

        f_nSRR = self.counters[key2.format(structures[1], neighbor, amino_acid, position)] + \
            self.counters[key2.format(
                structures[0], neighbor, amino_acid, position)]

        f_nSR = self.counters[key1.format(
            structures[0], amino_acid)] + self.counters[key1.format(structures[1], amino_acid)]

        return log10(f_SRR/f_nSRR) + log10(f_nSR/f_SR)

    def update_counters(self, sequence, structure):
        """
        Update the counters as asked : structure in which is the central amino acid, as well as its neighbourhood.
        :param sequence: str, amino acid sequence
        :param structure: str, sequence of the related structure
        :return: None, update the counters
        """
        neighborhood = [i for i in range(-8, 9) if i != 0]

        for i in range(len(sequence)):
            aa_in_struct = "{},{}".format(structure[i], sequence[i])

            if self.counters.get(structure[i]) == None:
                self.counters[structure[i]] = 0

            if self.counters.get(aa_in_struct) == None:
                self.counters[aa_in_struct] = 0

            self.counters[structure[i]] += 1
            self.counters[aa_in_struct] += 1

            for j in neighborhood:
                if i + j < 0:
                    continue
                elif i + j == len(sequence):
                    break

                neighbor_aa = "{},{},{},{}".format(
                    structure[i], sequence[i+j], sequence[i], str(j))

                if self.counters.get(neighbor_aa) == None:
                    self.counters[neighbor_aa] = 0

                self.counters[neighbor_aa] += 1

    def predict(self, sequence):
        """
        Compute the probability of the three structures for each amino acid according to the counters and
        select the highest as the prediction (I(S,R) + somme(I(S,R,Rj) with -8<=j<=8 and j=!0))
        :param sequence: str, sequence to predict
        :return: str, sequence of the predicted structure
        """

        structures = ["C", "H", "E"]
        neighborhood = [i for i in range(-8, 9) if i != 0]

        predicted_struct = ""
        for i in range(len(sequence)):
            probabilities = [0, 0, 0]
            for s in range(3):
                n_sum = 0

                for m in neighborhood:
                    if i + m < 0:
                        continue
                    elif i + m == len(sequence):
                        break

                    n_sum += self.neighborhood_inf(
                        sequence[i], structures[s], sequence[i+m], m)

                probabilities[s] = self.individual_inf(
                    sequence[i], structures[s]) + n_sum

            predicted_struct += structures[probabilities.index(
                max(probabilities))]

        return predicted_struct

    def validate(self, test_set):
        """
        Specific function to test with a test set of variable length
        :param test_set: list of tuples of size 2 (str, str), respectively amino acid sequence and structure sequence
        :return: tuple of size 4 each of size 2((float, float), (float, float), (float, float), (float, float)),
        respectively mean Q3, MCC-H, MMC-E and MCC-C and their standard deviation, all rounded at 2 decimal places
        """

        q_3s = np.array([], float)
        mcc_Hs = np.array([], float)
        mcc_Es = np.array([], float)
        mcc_Cs = np.array([], float)

        for sequence in test_set:
            predicted = self.predict(sequence[0])
            q_3s = np.append(q_3s, self.Q3(sequence[1], predicted))
            mcc_Hs = np.append(mcc_Hs, self.MCC(sequence[1], predicted, "H"))
            mcc_Es = np.append(mcc_Hs, self.MCC(sequence[1], predicted, "E"))
            mcc_Cs = np.append(mcc_Hs, self.MCC(sequence[1], predicted, "C"))

        def mean(array):
            return round(np.mean(array), 2)

        def std(array):
            return round(np.std(array), 2)

        return ((mean(q_3s), std(q_3s)), (mean(mcc_Hs), std(mcc_Hs)), (mean(mcc_Es), std(mcc_Es)), (mean(mcc_Cs), std(mcc_Cs)))

    def Q3(self, real, predicted):
        """
        Q3 computation
        :param real: str, real structure sequence
        :param predicted: str, predicted structure sequence
        :return: float rounded at 2 decimal, percentage of good predictions overall
        """
        correct = 0
        total = len(real)
        for i in range(total):
            if real[i] == predicted[i]:
                correct += 1

        return round(correct/total * 100, 2)

    def MCC(self, real, predicted, structure):
        """
        Matthew coefficient correlation. Evaluates the performance of a classifier (here, we use the binary variation)
        Needs to have the True Positives, True Negatives, False Positives and False Negatives computed according to
        the structure given in parameter. Returns None in case of the denominator is non valid (division by 0)
        :param real: str, real structure sequence
        :param predicted: str, predicted structure sequence
        :param structure: str, structure you want to evaluate the MCC of
        :return: float rounded at 2 decimal, MCC value, or None if division error
        """
        tp = 0
        tn = 0
        fp = 0
        fn = 0

        for i in range(len(real)):
            if predicted[i] == structure and real[i] == structure:
                tp += 1
            elif predicted[i] == structure and real[i] != structure:
                fp += 1
            elif predicted[i] != structure and real[i] != structure:
                tn += 1
            else:
                fn += 1

        denominator = (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
        if denominator == 0:
            return 0

        result = ((tp*tn)-(fp*fn))/sqrt(denominator) * 100
        return round(result, 2)
