"""Plot three different visualizations of the computed segment length distributions:

 - A single bar plot showing the expected distributions for one codomain size
 - A heatmap showing how the distributions behave for increasing hash function codomains
 - A line plot comparing three different codomain sizes

"""
import matplotlib.pyplot as plt
# plt.switch_backend('agg')  # Use this for SSH sessions
import seaborn as sns
import pandas as pd
import os
from matplotlib.colors import LogNorm


def parse_filename(filename):
    """ Extract parameters from filename with the pattern
    prediction_w=100_k=200_m=4096.txt
    """
    w, k, C = None, None, None

    for token in os.path.splitext(os.path.basename(filename))[0].split("_"):
        if "=" in token:
            param, value = token.split("=")
            if param == "w":
                w = int(value)
            elif param == "k":
                k = int(value)
            elif param == "m":
                C = int(value)
            else:
                print(f"Invalid token! cannot parse token {token}")
    return w, k, C


def df_from_files(data_files):
    all_dfs = []
    for data_file in data_files:
        w, k, C = parse_filename(data_file)
        df = pd.read_csv(data_file, names=["Probability"])
        df["Segment Length"] = range(1, len(df["Probability"]) + 1)
        df["w"] = w
        df["k"] = k
        df["C"] = C
        all_dfs.append(df)
    return pd.concat(all_dfs)


def visualize_distributions():
    df = df_from_files(snakemake.input.distribution_files)
    pivoted = pd.pivot_table(
        df,
        values='Probability',
        index=['C'],
        columns="Segment Length",
    )

    sns.set(font="DejaVu Sans")
    fig, ax = plt.subplots(figsize=(14, 10))
    # Use a logarithmic color scale, very loosely following:
    # https://stackoverflow.com/questions/36898008/seaborn-heatmap-with-logarithmic-scale-colorbar
    sns.heatmap(
        pivoted,
        linewidths=0.0,  # prevent lines in pdf export
        rasterized=True,  # prevent lines in pdf export
        norm=LogNorm(vmin=pivoted.min().min(), vmax=pivoted.max().max()),
        ax=ax,
    )
    ax.set_ylabel("C   ", rotation=0)

    ax.set_title("Segment length probabilities for w=50 and K=150")
    sns.despine()
    plt.savefig(snakemake.output.heatmap_path, bbox_inches='tight')
    plt.clf()

    sns.set(font="DejaVu Sans", style=sns.axes_style("whitegrid", {'grid.linestyle': '--'}))
    scale = 2.6
    fig, ax = plt.subplots(figsize=(scale*4, scale*2))

    selected = [1, 25, 50, 100]
    dfs = []
    for selector in selected:
        dfs.append(df[df["Segment Length"] == selector])

    filtered_df = pd.concat(dfs)
    filtered_df["Segment Length"] = pd.Categorical(filtered_df["Segment Length"])
    filtered_df["X"] = "Measurement"

    g = sns.lineplot(
        data=filtered_df,
        x="C",
        y="Probability",
        hue="Segment Length",
        style="X",
        markers=True,
        dashes=False,
        ax=ax,
    )
    for t in ax.get_legend().texts:
        if t.get_text() == "X":
            t.set_text("")

    sns.despine()

    plt.savefig(snakemake.output.lineplot_path, bbox_inches='tight')

    plt.clf()
    sns.set(
        font="DejaVu Sans",
        style=sns.axes_style("whitegrid", {'grid.linestyle': '--'}),
    )
    # scale = 3
    data = df[df["C"] == max(df["C"])]
    g = sns.FacetGrid(data, row="C", height=5, aspect=2)
    g = g.map(plt.bar, "Segment Length", "Probability", width=1.0,)
    plt.yscale("log")
    sns.despine()
    plt.savefig(snakemake.output.single_path, bbox_inches='tight')


if __name__ == "__main__":
    visualize_distributions()
