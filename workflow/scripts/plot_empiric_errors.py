"""Visualize the predictions for several codomain sizes vs. one specific empiric dataset.

This is used for the overview graphics showing multiple multiple predictions in shades of purple.
"""
import seaborn as sns
import pandas as pd
import os
import matplotlib.pyplot as plt
plt.switch_backend('agg')  # Use this for SSH sessions


def parse_filename(filename):
    """ Extract parameters from filename
    """
    w, k, q, C, canon, hf, gc, genome = None, None, None, None, None, None, None, None

    for token in os.path.splitext(os.path.basename(filename))[0].split("_"):
        if "=" in token:
            param, value = token.split("=")
            if param == "w":
                w = int(value)
            elif param == "q":
                q = int(value)
            elif param == "m":
                C = int(value)
            elif param == "k":
                k = int(value)
            elif param == "canon":
                canon = value
            elif param == "hf":
                hf = value
            elif param == "gc":
                gc = float(value)
            elif param == "genome":
                genome = value
            else:
                print(f"Invalid token! cannot parse token {token}")
    return w, k, C, q, canon, hf, gc, genome


def df_from_files(data_files):
    all_dfs = []
    for data_file in data_files:
        w, k, C, q, canon, hf, gc, genome = parse_filename(data_file)
        df = pd.read_csv(data_file, names=["Probability"])
        df["Segment Length"] = range(1, len(df["Probability"]) + 1)
        df["w"] = w
        df["k"] = k
        df["C"] = C
        all_dfs.append(df)
    return pd.concat(all_dfs)


def empiric_df_from_file(data_file):
    w, k, C, q, canon, hf, gc, genome = parse_filename(data_file)
    df = pd.read_csv(
        data_file,
        sep=" ",
        names=["Segment Length", "Count"],
        dtype={"Count": int},
    )
    df["w"] = w
    df["q"] = q
    df["hf"] = hf
    df["canon"] = canon
    df["Probability"] = df["Count"] / sum(df["Count"])
    return df


def plot_successive_errors():
    prediction_df = df_from_files(snakemake.input.small_distribution_files + snakemake.input.large_distribution_files)
    empiric_df = empiric_df_from_file(snakemake.input.empiric)

    scale = 3
    prediction_df["X"] = "Measurement"
    sns.set(
        font="DejaVu Sans",
        style=sns.axes_style("whitegrid", {'grid.linestyle': '--'}),
    )
    fig, ax = plt.subplots(figsize=(scale*4, scale*2))
    ax.bar(empiric_df["Segment Length"], empiric_df["Probability"], color="gray")

    custom_palette = sns.color_palette("BuPu", len(set(prediction_df["C"])))
    sns.set_palette(custom_palette)
    sns.lineplot(
        data=prediction_df,
        x="Segment Length",
        y="Probability",
        hue="C",
        style="X",
        markers=True,
        dashes=False,
        ax=ax,
        legend=False,
        alpha=0.7,
    )
    for (C, C_df) in prediction_df.groupby("C", as_index=False):
        ax.text(
            C_df["Segment Length"].iloc[-1] + 1,
            C_df["Probability"].iloc[-1],
            f"C = {C}",
            verticalalignment="center" if C != 30000 else "top",
            horizontalalignment="left",
            fontsize=7
        )
    sns.despine()
    print(empiric_df[empiric_df["Segment Length"] >= 100])
    plt.yscale("log")
    plt.savefig(snakemake.output.errors_path, bbox_inches='tight')


if __name__ == "__main__":
    plot_successive_errors()
