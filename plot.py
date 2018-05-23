# Makes a scatterplot from a file with an x-y pair on each line
import matplotlib.pyplot as plt
fn = 'out.txt.signal';
with open(fn) as f:
    lines = f.readlines()
    xs = []
    ys = []
    for line in lines:
        xs.append(int(line.split()[0]))
        ys.append(int(line.split()[1]))
    plt.scatter(xs, ys)
    plt.title('SV Detection with Unique Kmers')
    plt.xlabel('Required coverage')
    plt.ylabel('Number of SVs detected')
    plt.savefig('unique.png')
