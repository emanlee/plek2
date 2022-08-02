
import numpy as np
import pandas

from functions import filter_fasta
from functions import get_kmer
from functions import get_ORF
from functions import contact
from functions import prediction
from functions import init_data
from functions import output_acc


import  argparse
parser = argparse.ArgumentParser(description='Instruction of PLEK2(more details in README.txt):')
parser.add_argument('-i','--input_file', help='Input file in fasta format',required=True)
parser.add_argument('-m','--model',help='ve: vertebrate , pl = plant ', required=True)
# parser.add_argument('-o','--output_file',help='Output file in fasta format', required=True)
args = parser.parse_args()

# read the file

file = filter_fasta(args.input_file)
get_kmer('kmer_seqs')
get_ORF('seq_lines.fasta')
contact('kmer_6.txt','orf_length.txt')

print("loading data ...")

nums = init_data("features.txt")

##准备数据进行预测
dataframe = pandas.DataFrame(nums)
dataset = dataframe.values

x = dataset.astype(float)

X = np.zeros(shape=(len(x),1,5461,1))

for i in range(len(x)):
	tp = np.asarray(x[i])
	tp = np.resize(tp,(1,5461,1))
	X[i] = tp


print("Loading Model and predicting ...")

md = args.model

out = prediction(X,md)

print(out)

prediction = output_acc(out)

# out_file = open(args.output_file,'w+')







