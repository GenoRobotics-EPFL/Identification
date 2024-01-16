import matplotlib.pyplot as plt

#creating plots for the number of clusters
plt.figure()
plt.title("number of clusters vs coverage")
plt.xlabel("number of clusters retained")
plt.ylabel("Coverage [%]")
plt.plot([1, 3, 5, 7], [65.23, 87.75, 88.16, 88.7], '-o', label="rbcL")
plt.plot([1, 3, 5, 7], [80.46, 73.3, 93.75, 93.75], '-o', label="matK")
plt.plot([1, 3, 5, 7], [18.79, 47.36, 47.36, 47.36], '-o', label="psbA-trnH")
plt.plot([1, 3, 5, 7], [87.02, 87.05, 87.9, 89.1], '-o', label="ITS2")
plt.legend()
plt.savefig("cluster_number_coverage.png")

#rbcL : 126.1s, matK : 109.28s, psbA-trnH : 110.62s, ITS2 : 15032s
plt.figure()
plt.title("Number of clusters vs time")
plt.xlabel("Number of clusters retained")
plt.ylabel("Computation Time relative to a standard execution of 3 clusters [%]")
plt.plot([1, 3, 5, 7], [76.95, 100, 128.56, 151.27], '-o', label="rbcL")
plt.plot([1, 3, 5, 7], [64.29, 100, 122.03, 154.61], '-o', label="matK")
plt.plot([1, 3, 5, 7], [59.3, 100, 136.98, 153.79], '-o', label="psbA-trnH")
plt.plot([1, 3, 5, 7], [12.82, 100, 118.8, 133.9], '-o', label="ITS2")
plt.legend()
plt.savefig("cluster_number_time.png")


#creating plots for the length of k-mers
plt.figure()
plt.title("length of k-mers vs coverage")
plt.xlabel("length of k-mers")
plt.ylabel("Coverage [%]")
plt.plot([4, 8, 15, 30, 60], [82.00, 88.85, 87.75, 88.16, 83.09], '-o', label="rbcL")
plt.plot([4, 8, 15, 30, 60], [84.68, 71.75, 73.30, 79.18, 79.18], '-o', label="matK")
plt.plot([4, 8, 15, 30, 60], [57.43, 57.14, 47.36, 45.16, 45.16], '-o', label="psbA-trnH")
plt.plot([4, 8, 15, 30, 60], [87.90, 87.90, 87.05, 87.90, 87.90], '-o', label="ITS2")
plt.legend()
plt.savefig("length_k_mer_coverage.png")

#rbcL : 126.1s, matK : 109.28s, psbA-trnH : 110.62s, ITS2 : 15032s
plt.figure()
plt.title("length of k-mers vs time")
plt.xlabel("length of k-mers")
plt.ylabel("Computation Time relative to a standard execution of 15-mers [%]")
plt.plot([4, 8, 15, 30, 60], [356.93, 109.75, 100, 86.03, 85.31], '-o', label="rbcL")
plt.plot([4, 8, 15, 30, 60], [130.44, 116.73, 100, 75.57, 77.96], '-o', label="matK")
plt.plot([4, 8, 15, 30, 60], [308.80, 213.73, 100, 65.87, 69.65], '-o', label="psbA-trnH")
plt.plot([4, 8, 15, 30, 60], [102.24, 110.94, 100, 70.02, 68.12], '-o', label="ITS2")
plt.legend()
plt.savefig("length_k_mer_time.png")



#creating plots for the cluster distance
plt.figure()
plt.title("maximal distance within clusters vs coverage")
plt.xlabel("maximal distance within clusters")
plt.ylabel("Coverage [%]")
plt.plot([0.001, 0.01, 0.1, 0.3, 0.5], [83.23, 83.23, 87.75, 87.06, 87.89], '-o', label="rbcL")
plt.plot([0.001, 0.01, 0.1, 0.3, 0.5], [81.09, 81.09, 73.30, 85.08, 88.13], '-o', label="matK")
plt.plot([0.001, 0.01, 0.1, 0.3, 0.5], [30.21, 23.40, 47.36, 57.91, 57.91], '-o', label="psbA-trnH")
plt.plot([0.001, 0.01, 0.1, 0.3, 0.5], [88.07, 88.07, 87.90, 87.90, 87.90], '-o', label="ITS2")
plt.legend()
plt.savefig("cluster_distance_coverage.png")

#rbcL : 126.1s, matK : 109.28s, psbA-trnH : 110.62s, ITS2 : 15032s
plt.figure()
plt.title("maximal distance within clusters vs time")
plt.xlabel("maximal distance within clusters")
plt.ylabel("Computation Time relative to a standard execution of 0.1 distance allowed [%]")
plt.plot([0.001, 0.01, 0.1, 0.3, 0.5], [84.92, 84.16, 100, 168.95, 342.83], '-o', label="rbcL")
plt.plot([0.001, 0.01, 0.1, 0.3, 0.5], [83.70, 79.37, 100, 145.36, 126.85], '-o', label="matK")
plt.plot([0.001, 0.01, 0.1, 0.3, 0.5], [71.74, 69.81, 100, 330.18, 318.39], '-o', label="psbA-trnH")
plt.plot([0.001, 0.01, 0.1, 0.3, 0.5], [67.95, 70.01, 100, 68.52, 237.25], '-o', label="ITS2")
plt.legend()
plt.savefig("cluster_distance_time.png")


#creating plots for the mash distance
plt.figure()
plt.title("mash distance both sides vs coverage")
plt.xlabel("mash distance both sides")
plt.ylabel("Coverage [%]")
plt.plot([50, 100, 200, 400], [87.34, 87.27, 87.75, 85.49], '-o', label="rbcL")
plt.plot([50, 100, 200, 400], [92.73, 83.78, 73.30, 92.73], '-o', label="matK")
plt.plot([50, 100, 200, 400], [57.91, 49.76, 47.36, 40.65], '-o', label="psbA-trnH")
plt.plot([50, 100, 200, 400], [87.05, 87.05, 87.05, 87.90], '-o', label="ITS2")
plt.legend()
plt.savefig("mash_distance_coverage.png")

#rbcL : 126.1s, matK : 109.28s, psbA-trnH : 110.62s, ITS2 : 15032s
plt.figure()
plt.title("mash distance both sides vs time")
plt.xlabel("mash distance both sides")
plt.ylabel("Computation Time relative to a standard execution for 200 mash distance [%]")
plt.plot([50, 100, 200, 400], [106.25, 97.95, 100, 121.30], '-o', label="rbcL")
plt.plot([50, 100, 200, 400], [94.15, 95.42, 100, 124.47], '-o', label="matK")
plt.plot([50, 100, 200, 400], [274.86, 177, 100, 111.59], '-o', label="psbA-trnH")
plt.plot([50, 100, 200, 400], [24.17, 100, 100, 103.78], '-o', label="ITS2")
plt.legend()
plt.savefig("mash_distance_time.png")
