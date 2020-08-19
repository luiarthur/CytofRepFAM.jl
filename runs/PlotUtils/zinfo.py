import numpy as np
import matplotlib.pyplot as plt

def plot_num_selected_features(R, refs=None, ymax=1, 
                               xlabel="# of selected features"):
                               
    I, K = R.shape
    prop_selected = R / R.sum(1, keepdims=True)
    for i in range(I):
        plt.subplot(I, 1, i + 1)
        plt.plot(np.arange(K) + 1, prop_selected[i, :], marker=".", c="black")
        plt.xticks(np.arange(K) + 1, rotation=90)
        plt.ylabel(f"probability (sample {i + 1})")
        if refs is not None:
            plt.axvline(refs[i], ls=":", color="red")
        plt.ylim(-0.05, ymax)

    plt.xlabel(xlabel)


# TEST
# np.random.seed(1)
# R = np.random.randint(2, size=(2, 15))
# 
# plt.figure(figsize=(5,5))
# plot_num_selected_features(R, refs=[6, 5])
# plt.tight_layout()
# plt.show()
