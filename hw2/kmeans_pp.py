#import mykmeanssp as kmeans
import numpy as np
import pandas as pd
#import sys

def Kmeans_pp(datapoints,k):
    n = len(datapoints)
    centroids = [0]
    np.random.seed(0)
    centroids[0] = np.random.choice(n,1)[0]

    for z in range(k):
        # maybe can empty the existing array...?
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
    return data.drop(0,1)


def arrToSeq(arr):
    return [item for vec in arr for item in vec]


def main():
    k=3
    datapoints_table = readAndJoin('test_data\input_1_db_1.txt', 'test_data\input_1_db_2.txt')
    datapoints_metrix = datapoints_table.values.tolist()
    centroidArray = Kmeans_pp(datapoints_metrix, k)
    print(centroidArray)


# def main():
#     try:
#         k = int(sys.argv[1])
#         if (k <= 0):
#             print("INPUT ERROR:\nk can't be <= 0")
#             return 1
#     except ValueError:
#         print("INPUT ERROR:\nk can't be a letter")
#         return 1
#     if len(sys.argv) >= 5:
#         try:
#             max_iter = int(sys.argv[2])
#
#         except ValueError:
#             print("INPUT ERROR:\nk can't be a letter")
#             return 1
#
#         if max_iter <= 0:
#             print("INPUT ERROR:\nmax iteration can't be <= 0")
#             return 1
#         datapoints_table = readAndJoin(sys.argv[3], sys.argv[4])
#         datapoints_metrix = datapoints_table.values.tolist()
#         centroidArray = Kmeans_pp(datapoints_metrix, k)
#         datapoints = arrToSeq(datapoints_metrix)
#         centroid = arrToSeq(centroidArray.tolist())
#         centroid = kmeans.fit(k, datapoints, centroid, max_iter, len(datapoints_table[0]), len(datapoints_table))


if __name__ == "__main__":
    main()



