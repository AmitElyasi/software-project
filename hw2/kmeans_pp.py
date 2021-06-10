import numpy as np
import pandas as pd

def Kmeans_pp(datapoints,k):
    n = len(datapoints)
    centroids = np.empty(k)
    np.random.seed(0)
    centroids[0] = np.random.choice(range(n))

    for z in range(k):
        # maybe can empty the existing array...?
        min_dists = np.empty(n)

        for point in datapoints:
            print(type(datapoints))
            np.append(min_dists, min(np.inner(point-centroid) for centroid in centroids))

        dists_sum = np.sum(min_dists)
        probabilities = np.array(dist/dists_sum for dist in min_dists)
        new_centroid = np.random.choice(datapoints, p = probabilities)
        centroids.append(new_centroid[0])

    return centroids

def read_data(file1, file2=None):
    data1 = pd.read_csv(file1)
    if file2:
        data2 = pd.read_csv(file2)
    return data1



if __name__ == '__main__':
    datapoints = read_data('input_1_db_1.txt')
    #print(datapoints)
    print(Kmeans_pp(datapoints, 3))
