import os
import numpy as np
from tensorflow.keras.models import load_model

from itertools import chain
from Bio import Seq

import regex as re


################  filter_fasta file ##############################
def filter_fasta(fa):
    # Change the format of the input fasta file ...
    f1 = open(fa, 'r').readlines()
    f2 = open('seq_to_one_line.fasta', 'w')
    n = 0
    for i in f1:
        if i.startswith('>'):
            n += 1
            if n == 1:
                f2.write(i)
            else:
                f2.write('\n' + i)
        else:
            f2.write(i.strip("\n"))  # strip function strips the left and right Spaces
    f2.close()

    # Remove seqences whose lengths are not more than minlength ...

    file_in = open("seq_to_one_line.fasta", 'r')

    fa_Con = file_in.read()

    file_in.close()

    every_fas = fa_Con.split(">")

    f3 = open("filter_sequence_minlength.fasta", 'w')

    minlength = 200

    for i in every_fas:
        if i != "":
            start = i.index("\n")
            # print(i[start:])
            if len(i[start:]) > minlength:
                f3.write(">" + i)
    f3.close()

    # Split define lines and nucleotide lines

    file_in = open("filter_sequence_minlength.fasta", 'r').readlines()

    f4 = open("define_lines.fasta", 'w')
    f5 = open("seq_lines.fasta", 'w')

    n = 0
    for i in file_in:
        if i.startswith('>'):
            n += 1
            if n == 1:
                f4.write(i)
            else:
                f4.write('\n' + i)
        else:
            f5.write(i)

    f4.close()
    f5.close()

    # Nucleotides to upper case ...

    file_in = open("seq_lines.fasta", 'r').readlines()

    f6 = open("seq_upper.fasta", 'w')

    for i in file_in:
        f6.write(i.upper())

    f6.close()

    # Replace U with T

    file_in = open('seq_upper.fasta', 'r')
    f_new = open('replace_u', 'w')
    find_str = 'U'
    replace_str = 'T'
    for line in file_in:
        if find_str in line:
            line = line.replace(find_str, replace_str)
            f_new.write(line)
        else:
            f_new.write(line)

    f_new.close()

    # Replace R Y M K S W H B V D with N
    file_in = open('replace_u', 'r')
    f_new = open('replace_R', 'w')
    find_str = 'R'
    replace_str = 'N'
    for line in file_in:
        if find_str in line:
            line = line.replace(find_str, replace_str)
            f_new.write(line)
        else:
            f_new.write(line)

    f_new.close()

    file_in = open('replace_R', 'r')
    f_new = open('replace_Y', 'w')
    find_str = 'Y'
    for line in file_in:
        if find_str in line:
            line = line.replace(find_str, replace_str)
            f_new.write(line)
        else:
            f_new.write(line)

    f_new.close()

    file_in = open('replace_Y', 'r')
    f_new = open('replace_M', 'w')
    find_str = 'M'
    for line in file_in:
        if find_str in line:
            line = line.replace(find_str, replace_str)
            f_new.write(line)
        else:
            f_new.write(line)

    f_new.close()

    file_in = open('replace_M', 'r')
    f_new = open('replace_K', 'w')
    find_str = 'K'
    for line in file_in:
        if find_str in line:
            line = line.replace(find_str, replace_str)
            f_new.write(line)
        else:
            f_new.write(line)

    f_new.close()

    file_in = open('replace_K', 'r')
    f_new = open('replace_S', 'w')
    find_str = 'S'
    for line in file_in:
        if find_str in line:
            line = line.replace(find_str, replace_str)
            f_new.write(line)
        else:
            f_new.write(line)

    f_new.close()

    file_in = open('replace_S', 'r')
    f_new = open('replace_W', 'w')
    find_str = 'W'
    for line in file_in:
        if find_str in line:
            line = line.replace(find_str, replace_str)
            f_new.write(line)
        else:
            f_new.write(line)

    f_new.close()

    file_in = open('replace_W', 'r')
    f_new = open('replace_H', 'w')
    find_str = 'H'
    for line in file_in:
        if find_str in line:
            line = line.replace(find_str, replace_str)
            f_new.write(line)
        else:
            f_new.write(line)

    f_new.close()

    file_in = open('replace_H', 'r')
    f_new = open('replace_B', 'w')
    find_str = 'B'
    for line in file_in:
        if find_str in line:
            line = line.replace(find_str, replace_str)
            f_new.write(line)
        else:
            f_new.write(line)

    f_new.close()

    file_in = open('replace_B', 'r')
    f_new = open('replace_V', 'w')
    find_str = 'V'
    for line in file_in:
        if find_str in line:
            line = line.replace(find_str, replace_str)
            f_new.write(line)
        else:
            f_new.write(line)

    f_new.close()

    file_in = open('replace_V', 'r')
    f_new = open('kmer_seqs', 'w')
    find_str = 'D'
    for line in file_in:
        if find_str in line:
            line = line.replace(find_str, replace_str)
            f_new.write(line)
        else:
            f_new.write(line)

    f_new.close()

    ## remove middle file
    os.remove('replace_u')
    os.remove('replace_R')
    os.remove('replace_Y')
    os.remove('replace_M')
    os.remove('replace_K')
    os.remove('replace_S')
    os.remove('replace_W')
    os.remove('replace_H')
    os.remove('replace_B')
    # os.remove('replace_V')   此时还不能删除replace_V


#################### get kmer_features #####################################
def get_kmer(kmer_seqs):
    class kmer_featurization:

        def __init__(self, k):
            """
            seqs: a list of DNA sequences
            k: the "k" in k-mer
            """
            self.k = k
            self.letters = ['A', 'T', 'C', 'G']
            self.multiplyBy = 4 ** np.arange(k - 1, -1,
                                             -1)  # the multiplying number for each digit position in the k-number system
            self.n = 4 ** k  # number of possible k-mers

        def obtain_kmer_feature_for_a_list_of_sequences(self, seqs, write_number_of_occurrences=False):
            """
            Given a list of m DNA sequences, return a 2-d array with shape (m, 4**k) for the 1-hot representation of the kmer features.

            Args:
              write_number_of_occurrences:
                a boolean. If False, then in the 1-hot representation, the percentage of the occurrence of a kmer will be recorded; otherwise the number of occurrences will be recorded. Default False.
            """
            kmer_features = []
            for seq in seqs:
                this_kmer_feature = self.obtain_kmer_feature_for_one_sequence(seq.upper(),
                write_number_of_occurrences=write_number_of_occurrences)
                kmer_features.append(this_kmer_feature)

            kmer_features = np.array(kmer_features)

            return kmer_features

        def obtain_kmer_feature_for_one_sequence(self, seq, write_number_of_occurrences=False):
            """
            Given a DNA sequence, return the 1-hot representation of its kmer feature.

            Args:
              seq:
                a string, a DNA sequence
              write_number_of_occurrences:
                a boolean. If False, then in the 1-hot representation, the percentage of the occurrence of a kmer will be recorded; otherwise the number of occurrences will be recorded. Default False.
            """
            number_of_kmers = len(seq) - self.k + 1
            kmer_feature = np.zeros(self.n)

            # k-mer features is transformed into one-hot encoding
            for i in range(number_of_kmers):
                temporary = seq[i:(i + self.k)]
                if 'N' not in temporary:  # if 'N' in kmer,no action is performed.(skip 'N')
                    this_kmer = seq[i:(i + self.k)]
                    this_numbering = self.kmer_numbering_for_one_kmer(this_kmer)
                    kmer_feature[this_numbering] += 1

            if not write_number_of_occurrences:
                kmer_feature = (kmer_feature / number_of_kmers) * pow(4, (self.k) - 5)

            return kmer_feature

        def kmer_numbering_for_one_kmer(self, kmer):
            """
            Given a k-mer, return its numbering (the 0-based position in 1-hot representation)
            display the position of k
            """
            digits = []
            for letter in kmer:
                digits.append(self.letters.index(letter))
            digits = np.array(digits)
            numbering = (digits * self.multiplyBy).sum()
            return numbering

    ### get k=6 features
    # Assign the file to the list line by line, removing the newline after each line
    with open(kmer_seqs, 'r') as f:
        seq_list = f.read().splitlines()

    # file line count, count how many samples
    count_line = 0
    f = open(kmer_seqs, 'r')
    for line in f.readlines():
        count_line = count_line + 1

    k = 1
    obj = kmer_featurization(k)
    kmer_feature_1 = obj.obtain_kmer_feature_for_a_list_of_sequences(seq_list, write_number_of_occurrences=False)
    kmer_feature_1 = kmer_feature_1.reshape(count_line, 4)

    k = 2
    obj = kmer_featurization(k)
    kmer_feature_2 = obj.obtain_kmer_feature_for_a_list_of_sequences(seq_list, write_number_of_occurrences=False)
    kmer_feature_2 = kmer_feature_2.reshape(count_line, 16)

    k = 3
    obj = kmer_featurization(k)
    kmer_feature_3 = obj.obtain_kmer_feature_for_a_list_of_sequences(seq_list, write_number_of_occurrences=False)
    kmer_feature_3 = kmer_feature_3.reshape(count_line, 64)

    k = 4
    obj = kmer_featurization(k)
    kmer_feature_4 = obj.obtain_kmer_feature_for_a_list_of_sequences(seq_list, write_number_of_occurrences=False)
    kmer_feature_4 = kmer_feature_4.reshape(count_line, 256)

    k = 5
    obj = kmer_featurization(k)
    kmer_feature_5 = obj.obtain_kmer_feature_for_a_list_of_sequences(seq_list, write_number_of_occurrences=False)
    kmer_feature_5 = kmer_feature_5.reshape(count_line, 1024)

    k = 6
    obj = kmer_featurization(k)
    kmer_feature_6 = obj.obtain_kmer_feature_for_a_list_of_sequences(seq_list, write_number_of_occurrences=False)
    kmer_feature_6 = kmer_feature_6.reshape(count_line, 4096)

    # k=6
    kmer_6 = np.concatenate(
        (kmer_feature_1, kmer_feature_2, kmer_feature_3, kmer_feature_4, kmer_feature_5, kmer_feature_6), axis=1)
    np.savetxt("kmer_6.txt", kmer_6)

def get_ORF(fa):
    # read data
    def read_data(seq_line):
        RNA_data = []
        try:
            with open(seq_line, "rt") as fp:
                lines = fp.readlines()
            for i in range(0, len(lines)):
                RNA_data.append(lines[i].replace("\n", "").strip().split())
        except:
            print("抛出异常")
        finally:
            return RNA_data

    seq = read_data(fa)
    seq = list(chain.from_iterable(seq))

    ###   Calculate the length of the peptide chain,length is the length of the peptide chain
    startP = re.compile('ATG')  # The starting position is ATG
    nuc = seq
    longest = (0,)  # The length of the peptide chain
    length = []
    for i in range(len(nuc)):
        for m in startP.finditer(nuc[i], overlapped=True):  # find the location of all the start codons
            if len(Seq.Seq(nuc[i])[m.start():].translate(to_stop=True)) > longest[0]:  # longet[0] , the starting position of ORF is the longest length of the peptide chain
                pro = Seq.Seq(nuc[i])[m.start():].translate(to_stop=True)  # pro = peptides
                longest = (len(pro),
                           m.start(),
                           str(pro),
                           nuc[i][m.start():m.start() + len(pro) * 3 + 3])
        longest = (0,)  # at the end of each large cycle, the length of the peptide chain is assigned a value of 0
        length.append(len(pro) + 1)

    # Take the log base 10 of the length of the peptide chain
    length_log10 = []
    length_log10 = np.log10(length)

    orf_normalize = []
    for i in range(len(length_log10 - 1)):
        x = float(length_log10[i] - 0.48) / (4 - 0.48) * 0.1
        orf_normalize.append(x)

    ### Save the list data in TXT format by row and read it
    with open('orf_length.txt', 'w') as f:
        for i in orf_normalize:
            f.write(str(i) + '\n')

def contact(kmer,ORF):
    k_mer = np.loadtxt(kmer)
    ORF = np.loadtxt(ORF)

    features = np.column_stack((k_mer, ORF))
    np.savetxt('features.txt', features, fmt='%0.8f')

    os.remove('replace_V')
    os.remove('seq_upper.fasta')

###################################### The feature file is input into the network model for prediction   ##################################

def read_data(human_data):
   RNA_data = []
   try:
      with open(human_data,"rt") as fp:
         lines = fp.readlines()
         # print(lines)
      for i in range(0,len(lines)):           # Read from the first line
         RNA_data.append(lines[i].replace("\n","").strip().split())
   except:
      print("抛出异常")
   finally:
      return RNA_data


## 初始化数据
def init_data(human_data):
   temp_datas = read_data(human_data)
   if len(temp_datas)<=0:
      print("数据初始化失败")
   else:
      RNA_data = np.array(temp_datas)     #create array
      nums = RNA_data.astype(np.float)
      return nums

nums = init_data("features.txt")

def prediction(dat,md):
   if md == 've':
       model = load_model('Coding_Net_kmer6_orf.h5')
   else:
       model = load_model('Coding_Net_kmer6_orf_Arabidopsis.h5')
   Y = model.predict_classes(dat)
   labels = np.array(Y)
   # save results file
   np.savetxt('results',labels,fmt='%d')

   return labels


def output_acc(Y):
   seq_len = len(Y)
   noncoding = 0
   coding = 1
   for i in range(len(Y)):
      if Y[i] == 0:
         noncoding = noncoding + 1
      else:
          coding = coding + 1
   pred_noncoding_acc = noncoding / seq_len
   pred_coding_acc = 1-pred_noncoding_acc
   acc = [pred_noncoding_acc,pred_coding_acc]
   print('non-coding =', pred_noncoding_acc)
   print('coding =' , pred_coding_acc)

   # remove middle file
   os.remove('define_lines.fasta')
   os.remove('filter_sequence_minlength.fasta')
   os.remove('kmer_6.txt')
   os.remove('kmer_seqs')
   # os.remove('orf_length.txt')
   os.remove('seq_lines.fasta')
   os.remove('seq_to_one_line.fasta')

   return acc



