# Makes a scatterplot from a file with an x-y pair on each line
import matplotlib.pyplot as plt
import sys

fn = sys.argv[1]
pr = False
with open(fn) as f:
    lines = f.readlines()
    xs = []
    ys = []
    for line in lines:
        xs.append(int(line.split()[0]))
        ys.append(int(line.split()[1]))
    plt.scatter(xs, ys)
    if not pr:
        plt.title('SV Detection with Unique Kmers')
        plt.xlabel('Required coverage')
        plt.ylabel('Number of SVs detected')
    else:
        plt.title('Precision/Recall for Unique Kmer SV Detection')
        plt.xlabel('Precision')
        plt.ylabel('Recall')
    plt.savefig(sys.argv[2])
