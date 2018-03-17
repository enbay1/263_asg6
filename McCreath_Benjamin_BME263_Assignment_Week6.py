import csv
import math
import sys
import time
from operator import itemgetter

import matplotlib.patches as mplpatches
import matplotlib.pyplot as plt
from numpy import arange


def make_circle(r=0, range_start=0, range_end=360):
    """Return the x and y points of a circle or arc with radius r."""
    if range_start > range_end:
        temp = range_start
        range_start = range_end
        range_end = temp
    return [[r * math.cos(angle * (math.pi / 180)),
             r * math.sin(angle * (math.pi / 180))] for angle in arange(range_start, range_end, .1)]


def load_data() -> list:
    """Open data file and read it in."""
    try:
        with open(sys.argv[1]) as data_file:
            reader = csv.reader(data_file, delimiter='\t')
            next(reader)
            data = list(reader)
        return data
    # Catches user error.
    except FileNotFoundError as e:
        if isinstance(e, FileNotFoundError):
            sys.stderr.write("One of the specified files was not found.\n")
        else:
            sys.stderr.write("Unexpected error\n")
        sys.exit(1)


def process_data(data: list) -> (list, list):
    """Process raw data into data for right and left plots."""
    data = [itemgetter(*list(range(4, 12)), 13)(row) for row in data]
    data = [list(map(float, row)) for row in data]
    data = sorted(data, key=itemgetter(-1))
    left_data = []
    right_data = []
    for line in data:
        left_data.append(line[:-1])
        right_data.append(line[-1])
    return left_data, right_data


def get_binned_data(data: list) -> dict:
    """Bin data based on value, to be plotted on circular histogram."""
    data = sorted(data)
    binned_data = {}
    j = 0
    for i in range(2, 25, 2):
        binned_data[i] = []
        while j < len(data) and data[j] < i:
            binned_data[i].append(data[j])
            j += 1
    for data_bin in binned_data:
        binned_data[data_bin] = len(binned_data[data_bin])
    return binned_data


def normalize_radial_data(data: int) -> float:
    """Normalize data into circular histagram data-space."""
    return (data * 225) / 300 + 75


def main():
    """Drive plot creation and set properties. Load data and such as necessary."""
    start = time.time()
    data = load_data()
    left_data, right_data = process_data(data)
    del data
    right_data = get_binned_data(right_data)
    plt.style.use('BME163.mplstyle')
    plt_size = [5, 3]
    plt.figure(figsize=plt_size)
    left_panel = plt.axes([0.5 / plt_size[0], 0.3 / plt_size[1], 0.75 / plt_size[0], 2.5 / plt_size[1]],
                          xlim=(-1.5, 22.5), ylim=(0, 1260), xlabel="CT", ylabel="Number of genes",
                          xticks=list(range(0, 22, 3)),
                          xticklabels=[str(intlab) if not intlab % 6 else "" for intlab in range(0, 22, 3)])
    right_panel = plt.axes([1.75 / plt_size[0], 0.297 / plt_size[1], 2.5 / plt_size[0], 2.5 / plt_size[1]],
                           frame_on=False, xlim=[-300, 300], ylim=[-300, 300], xticks=[], yticks=[])
    del plt_size
    ################################ Right Panel #######################################################################
    # Draw dashed grid lines
    for num in (150, 225, 300):
        circle = zip(*make_circle(num))
        right_panel.plot(next(circle), next(circle), lw=0.3, dashes=[4, 8, 8, 8],
                         color="black")
        right_panel.text(-1 * num, 0, str(int(num - 50 * (1 / (num // 100)))) if num is not 300 else str(300),
                         fontsize=6, horizontalalignment="right", verticalalignment="center")
    # Draw inner circle
    for num in range(60, 76):
        circle = zip(*make_circle(num, range_start=90 if num not in [60, 75] else 0,
                                  range_end=270 if num not in [60, 75] else 360))
        right_panel.plot(next(circle), next(circle), lw=0.3 if num in [60, 75] else .5, color="black")
    # Draw inner text
    for i, angle in enumerate(range(90, 430, 60)):
        radius = 37.5
        right_panel.text(-radius * math.cos(angle * (math.pi / 180)), radius * math.sin(angle * (math.pi / 180)),
                         str(i * 4), fontsize=6, horizontalalignment="center", verticalalignment="center")
    right_panel.text(0, 0, "CT", fontsize=6, horizontalalignment="center", verticalalignment="center")
    # Draw data blocks
    for data_bin in right_data:
        base_one = int((data_bin - 2) * (-360 / 24) + 90)
        base_two = int(data_bin * (-360 / 24) + 90)
        r = normalize_radial_data(right_data[data_bin])
        for i in range(76, int(r)):
            circle = zip(*make_circle(i, base_two, base_one))
            right_panel.plot(next(circle), next(circle), lw=.4, color=[128 / 255] * 3, zorder=0)
        circle = zip(*make_circle(r, base_two, base_one))
        right_panel.plot(next(circle), next(circle), lw=.5, color="black")
        # draw "vertical" lines
        lower = itemgetter(0, -1)(make_circle(75, base_two, base_one))
        upper = itemgetter(0, -1)(make_circle(r, base_two, base_one))
        for point in zip(lower, upper):
            right_panel.plot([x[0] for x in point], [y[1] for y in point], lw=.5, color="black")
    ################################# Left Panel #######################################################################
    for i, row in enumerate(reversed(left_data)):
        for j, value in enumerate(row):
            plot_value = ((value - min(row)) / (max(row) - min(row))) * 100
            max_blue = [56, 66, 157]
            max_yellow = [253, 223, 41]
            norm_plot_v = plot_value / 100
            color = [(b * norm_plot_v + y * (1 - norm_plot_v)) / 255 for b, y in zip(max_blue, max_yellow)]
            left_panel.add_patch(mplpatches.Rectangle([j * 3 - 1.5, i], 3, 1, facecolor=color, linewidth=0))
    plt.savefig('McCreath_Benjamin_BME263_Assignment_Week6.png', dpi=600)
    print("Time to complete: {}s".format(time.time() - start))


if __name__ == '__main__':
    main()
