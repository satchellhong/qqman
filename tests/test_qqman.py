import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from qqman import qqman  # noqa: E402


def make_assoc(chroms, per_chrom=200):
    rng = np.random.default_rng(0)
    n = per_chrom * len(chroms)
    return pd.DataFrame(
        {
            "CHR": np.repeat(chroms, per_chrom),
            "BP": np.tile(np.arange(per_chrom) * 100, len(chroms)),
            "P": rng.uniform(1e-9, 1, n),
            "SNP": [f"rs{i}" for i in range(n)],
        }
    )


@pytest.mark.parametrize(
    "chroms, labels",
    [([3, 1, 2], ["1", "2", "3"]), (["2", "X", "1"], ["1", "2", "X"])],
)
def test_manhattan_plots_every_snp_in_chromosome_order(chroms, labels):
    df = make_assoc(chroms)
    fig, ax = plt.subplots()

    qqman.manhattan(df, ax=ax)

    assert [t.get_text() for t in ax.get_xticklabels()] == labels
    plotted = np.concatenate([c.get_offsets()[:, 1] for c in ax.collections])
    expected = -np.log10(df.sort_values(["CHR", "BP"])["P"].to_numpy())
    np.testing.assert_allclose(plotted, expected)
    plt.close(fig)


def test_manhattan_reads_whitespace_delimited_file(tmp_path):
    path = tmp_path / "assoc.txt"
    make_assoc([1, 2]).to_csv(path, sep=" ", index=False)
    fig, ax = plt.subplots()

    qqman.manhattan(str(path), ax=ax)

    assert sum(len(c.get_offsets()) for c in ax.collections) == 400
    plt.close(fig)


@pytest.mark.parametrize("kind", ["dataframe", "series", "list"])
def test_qqplot_drops_out_of_range_p_values(kind):
    df = make_assoc([1, 2])
    df.loc[:4, "P"] = [0.0, 1.0, np.nan, 1.5, -0.1]
    assoc = {"dataframe": df, "series": df["P"], "list": list(df["P"])}[kind]
    fig, ax = plt.subplots()

    qqman.qqplot(assoc, ax=ax)

    observed = ax.collections[0].get_offsets()[:, 1]
    assert len(observed) == 395
    np.testing.assert_allclose(
        observed, -np.log10(np.sort(df["P"].iloc[5:].to_numpy()))
    )
    plt.close(fig)
