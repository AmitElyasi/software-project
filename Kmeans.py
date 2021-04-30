import numpy as np

def Kmeans(k,filename, max_iter=200):
    data_points = read_data(filename)
    if len(data_points <= k):
        print(f"there is less then k={k} data points")
        return

    centroids = data_points[:k]
    for i in range[max_iter]:
        clusters = assign(data_points, centroids)
        new_centroids = re_estimate(clusters, centroids)
        if set(new_centroids) == set(centroids):
            break
    return clusters


# def read_data(filename):
#     with open (filename, rb) as f:
#         for line in f.readlines():
#             data_points = np.array(line.rstrip().split(','))
#     return data_points


def assign(data_points, centroids):
    clusters = {}
    for point in data_points:
        closest_centroind = centroids[0]
        min = np.linalg.norm(point - closest_centroind)

        for centroid in centroids[0:]:
            dst = np.linalg.norm(point - centroid)
            if dst <= min:
                min = dst
                closest_centroind = centroid

        if clusters[str(closest_centroind)]:
            clusters[str(closest_centroind)].append(point)
        else:
            clusters[str(closest_centroind)] = [point]
    return clusters


def re_estimate(clusters, centroinds):
    centroinds = []
    for centroid in centroinds:
        for point in clusters[str(centroid)]:
            sum += point
        centroinds.append(sum/len(clusters[str(centroid)]))
    return centroids
