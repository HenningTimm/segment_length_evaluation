"""Visualize predicted versus empiric segment length distributions
for reference genomes and simulated genomes.

As opposed to the landscape plots, these plots all show the same
canonization strategy, one genome per row und a elected hash
function per column.
"""
from matplotlib.patches import Rectangle as Rect
import seaborn as sns
import pandas as pd
import os
import colorbrewer
import colorsys
import matplotlib.pyplot as plt
plt.switch_backend('agg')  # Use this for SSH sessions

LIN_COLOR = [channel / 255 for channel in colorbrewer.PuBuGn[5][2]]
DARK_LIN_COLOR = [(channel - 40) / 255 for channel in colorbrewer.PuBuGn[5][2]]
LOG_COLOR = [channel / 255 for channel in colorbrewer.PuBuGn[5][4]]
LIGHT_LOG_COLOR = [max((channel + 10), 0) / 255 for channel in colorbrewer.PuBuGn[5][4]]
DARK_LOG_COLOR = [max((channel - 10), 0) / 255 for channel in colorbrewer.PuBuGn[5][4]]


def parse_filename(filename):
    """ Extract parameters from filename
    """
    w, q, C, canon, hf, gc, genome, k = None, None, None, None, None, None, None, None

    for token in os.path.splitext(os.path.basename(filename))[0].split("_"):
        if "=" in token:
            param, value = token.split("=")
            if param == "w":
                w = int(value)
            elif param == "q":
                q = int(value)
            elif param == "m":
                C = int(value)
            elif param == "canon":
                canon = value
            elif param == "hf":
                hf = value
            elif param == "gc":
                gc = float(value)
            elif param == "genome":
                genome = value
            elif param == "k":
                genome = k
            else:
                print(f"Invalid token! cannot parse token {token}")
    return w, q, C, canon, hf, gc, genome, k


def prediction_df_from_file(data_file):
    w, _, C, _, _, _, _, k = parse_filename(data_file)
    df = pd.read_csv(data_file, names=["Probability"])
    df["Segment Length"] = range(1, len(df["Probability"]) + 1)
    df["w"] = w
    df["k"] = k
    df["C"] = C
    return df


def df_from_files(data_files):
    prediction = prediction_df_from_file(snakemake.input.predicted_distribution)
    # rename to avoid conflicts with other dfs
    prediction = prediction.rename(columns={
        "Probability": "Expected Probability"
    })
    # drop all unneeded columns
    prediction = (prediction[["Segment Length", "Expected Probability"]])

    all_dfs = []
    for data_file in data_files:
        w, q, C, canon, hf, gc, genome, k = parse_filename(data_file)
        df = pd.read_csv(data_file, sep=" ", names=["Segment Length", "count"])
        # Add in expected distribution for each position.
        # This leaves NaN for all segment lengths >w which are replaced later
        df = pd.merge(df, prediction, on="Segment Length", how="outer")
        df["count"] = df["count"].fillna(0)
        df["w"] = w
        df["q"] = q
        df["canonicity"] = canon
        df["hf"] = hf
        df["genome"] = genome
        df["prob"] = df["count"]/sum(df["count"])

        all_dfs.append(df)
    return (pd.concat(all_dfs), w, q)


def sum_and_scatter(df, w):
    """Sum up all segments above w into one position (w+5)
    """

    beyond_K = pd.read_csv(snakemake.input.predicted_error, names=["Expected Probability"])

    # split dataframe at w
    above_w = df[df["Segment Length"] > w]
    w_and_lower = df[df["Segment Length"] <= w]
    summed = above_w.groupby(["hf", "canonicity", "w", "q", "genome"], as_index=False).aggregate(
        {
            "count": ["sum"],
            "prob": ["sum"],
            "Expected Probability": ["sum"],
        },
    ).reset_index()
    # Set dummy segment length to aggregate al values above w
    dummy_offset = 5
    summed["Segment Length"] = w + dummy_offset
    # Remove multiindex
    summed.columns = [col[0] for col in summed.columns]
    summed.reset_index()

    # Add probabilities for >K
    summed["Expected Probability"] = summed["Expected Probability"] + beyond_K["Expected Probability"][0]

    combined = pd.concat((w_and_lower, summed), sort=True)

    return combined, above_w


def plot_genome_comparison():

    df, w, q = df_from_files(snakemake.input.data_files)

    # replace NaNs introduced by the outer merge with the prediction with zeros
    # NaNs occur, when segment lengths above K occur
    df = df.fillna(0)

    # split values above w and below (including) w
    # for the values below w summ up all values above w
    # into on entry at w + 5
    summed, above_w = sum_and_scatter(df, w)

    summed["Expected Prob Line"] = summed["Expected Probability"].where(summed["Segment Length"] <= w)
    # plot empiric distributions
    heights = {
        30: 4,
        50: 5,
        100: 6,
    }
    sns.set(
        font="DejaVu Sans",
        style=sns.axes_style("whitegrid", {'grid.linestyle': '--'}),
        font_scale=2,
    )

    g = sns.FacetGrid(
        summed,
        row="genome",
        col="hf",
        height=heights[w],
        aspect=2,
        margin_titles=True,
    )
    g.fig.suptitle(f"Segment Length Distributions for $w={w}$, $q={q}$", y=1.01)
    colors = {
        "non": sns.color_palette()[0],
        "min": sns.color_palette()[1],
        "max": sns.color_palette()[2],
        }

    g1 = g.map(plt.bar, "Segment Length", "prob", color=colors[snakemake.wildcards.canonical])

    # fix broken z-order
    # save lower layer (layer 1, barplot)
    # so that later layers can be put above
    backgroundartists = []
    for ax in g1.axes.flat:
        for li in ax.lines + ax.collections:
            li.set_zorder(1)
            backgroundartists.append(li)
        beyond_w_bar = [rect for rect in ax.get_children() if isinstance(rect, Rect)][-2]
        # take the color of the bar andmake it darker
        (col_r, col_g, col_b, col_a) = beyond_w_bar.get_fc()
        col_h, col_l, col_s = colorsys.rgb_to_hls(col_r, col_g, col_b)
        col_r, col_g, col_b = colorsys.hls_to_rgb(col_h, col_l-0.2, col_s)
        beyond_w_bar.set_color((col_r, col_g, col_b, col_a))

    # plot predicted points and manually place them on a higher layer than the bars
    g2 = g.map(sns.scatterplot, "Segment Length", "Expected Probability", color="black", s=50)
    for ax in g2.axes.flat:
        for li in ax.lines + ax.collections:
            if li not in backgroundartists:
                li.set_zorder(5)
    g.map(
        sns.lineplot,
        "Segment Length",
        "Expected Prob Line",
        color="black",
        alpha=0.7,
        palette=sns.color_palette("Set2_r"),
    )
    # Adjust ticks so that the summed values at the dummy offset are labeled
    # correctly. currently the dummy value is at w + 5
    for ax in g2.axes.flat:
        labels = [item.get_text() if item.get_text()!=f"{w+5}" else ">w" for item in ax.get_xticklabels()]
        ax.set_xticklabels(labels)

        label = ax.get_ylabel()
        ax.set_ylabel(label if label != "Expected Prob Line" else "Probability")

    g.set(yscale="log")
    sns.despine()
    plt.savefig(snakemake.output.genome_comp_pdf, bbox_inches='tight')


if __name__ == "__main__":
    plot_genome_comparison()
