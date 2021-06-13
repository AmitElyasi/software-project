import mykmeanssp as kmeans
import numpy as np
import pandas as pd
import sys

def Kmeans_pp(datapoints,k):
    n = len(datapoints)
    centroids = [0]
    centroids[0] = np.random.choice(n,1)[0]

    for z in range(k-1):
        min_dists = np.empty(n)

        for point in datapoints:
            np.append(min_dists, min(np.inner(point-centroid,point-centroid) for centroid in centroids))

        dists_sum = np.sum(min_dists)
        probabilities = [dist/dists_sum for dist in min_dists]
        new_centroid_indx = np.random.choice(n, 1, p=probabilities)[0]
        centroids.append(new_centroid_indx)
    return centroids

def readAndJoin(file1, file2):
    data1 = pd.read_csv(file1, header=None)
    data2 = pd.read_csv(file2, header=None)
    data = pd.merge(data1, data2, on=0, sort=True)
    data.to_string('myData.txt')
    return data.drop(0,1)


def arrToSeq(arr):
    return [item for vec in arr for item in vec]


def seqRoundAndPrint(seq, d):
    n = len(seq)//d
    seq = np.round(seq, 4)
    for i in range(n):
        for j in range(d-1):
            print(seq[i*d+j],",",end="")
        print(seq[i * d + d-1])



# def main():
#     k = 3
#     datapoints_table = readAndJoin('test_data\input_1_db_1.txt', 'test_data\input_1_db_2.txt')
#     datapoints_matrix = datapoints_table.values.tolist()
#     centroidArray = Kmeans_pp(datapoints_matrix, k)
#     print(centroidArray)


def main():
    np.random.seed(0)

    try:
        k = int(sys.argv[1])
        if (k <= 0):
            print("INPUT ERROR:\nk can't be <= 0")
            return 1
    except ValueError:
        print("INPUT ERROR:\nk can't be a letter")
        return 1
    if len(sys.argv) >= 5:
        try:
            max_iter = int(sys.argv[2])

        except ValueError:
            print("INPUT ERROR:\nk can't be a letter")
            return 1

        if max_iter <= 0:
            print("INPUT ERROR:\nmax iteration can't be <= 0")
            return 1

        datapoints_table = readAndJoin(sys.argv[3], sys.argv[4])
        datapoints_matrix = datapoints_table.values.tolist()
        if len(datapoints_matrix) < k:
            print(f"INPUT ERROR:\nthere are less then k={k} data points")
            return False
        centroids_indexes = Kmeans_pp(datapoints_matrix, k)
        datapoints = arrToSeq(datapoints_matrix)
        centroids = arrToSeq([datapoints_matrix[i] for i in centroids_indexes])
        centroids = kmeans.fit(k, datapoints, centroids, max_iter, len(datapoints_matrix[0]), len(datapoints_matrix))
        print(",".join(str(indx) for indx in centroids_indexes))
        seqRoundAndPrint(centroids,len(datapoints_matrix[0]))

if __name__ == "__main__":
    main()



