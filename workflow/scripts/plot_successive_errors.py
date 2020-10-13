"""Visualize the descreasing errors for for increasing hf codomain sizes.
"""
import seaborn as sns
import pandas as pd
import os
import matplotlib.pyplot as plt
# plt.switch_backend('agg')  # Use this for SSH sessions


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
        df["Successive Errors"] = df["Probability"] - all_dfs[-1]["Probability"] if len(all_dfs) > 0 else 0
        df["Successive Difference Squares"] = df["Successive Errors"]**2
        all_dfs.append(df)
    return pd.concat(all_dfs)


def plot_successive_errors():
    df = df_from_files(snakemake.input.distribution_files)
    C_df = df.groupby("C", as_index=False).agg({"Successive Difference Squares": ["sum"]})

    # NOTE: the [1:] is used to trim of the leading zero
    # for the first value that has no predecessor to compute
    # a successive error
    Cs = C_df["C"][1:]
    squared_sums = C_df["Successive Difference Squares"]["sum"][1:]
    plotting_df = pd.DataFrame(
        {
            "C": Cs,
            "Successive Difference Squares": squared_sums,
            "X": "Measurement",
        }
    )

    scale = 3
    sns.set_style("whitegrid", {'grid.linestyle': '--'})
    fig, ax = plt.subplots(figsize=(scale*5, scale*2))


    sns.lineplot(
        data=plotting_df,
        x="C",
        y="Successive Difference Squares",
        style="X",
        markers=True,
        dashes=False,
        ax=ax,
    )
    for t in ax.get_legend().texts:
        if t.get_text() == "X":
            t.set_text("")

    plt.yscale("log")
    sns.despine()
    ax.set_title("Sum of squared differences between successive values for C")

    plt.savefig(snakemake.output.errors_path, bbox_inches='tight')


if __name__ == "__main__":
    plot_successive_errors()
