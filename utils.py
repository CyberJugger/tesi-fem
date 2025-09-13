import matplotlib.pyplot as plt

#prende come input le ascisse e la funzione
def plotter(xx, f, titolo_grafico, marker=False):
    plt.figure()
    if marker:
        plt.plot(xx, f, marker='o', markersize=1)
    else:
        plt.plot(xx, f)
    plt.xlabel("x"); plt.ylabel("T(x)")
    plt.title("{}".format(titolo_grafico))
    plt.grid(True)

