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


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def test_manhattan_leaves_input_dataframe_unchanged():
    df = make_assoc([1, 2])
    before = df.copy()
    fig, ax = plt.subplots()

    qqman.manhattan(df, ax=ax)

    pd.testing.assert_frame_equal(df, before)


@pytest.mark.parametrize("plot", [qqman.manhattan, qqman.qqplot])
def test_out_saves_the_figure_that_owns_ax(plot, tmp_path):
    fig, ax = plt.subplots()
    plt.subplots()  # a second, empty figure becomes pyplot's current figure

    plot(make_assoc([1, 2]), ax=ax, out=tmp_path / "plot.png")

    fig.savefig(tmp_path / "expected.png")
    np.testing.assert_array_equal(
        plt.imread(tmp_path / "plot.png"), plt.imread(tmp_path / "expected.png")
    )


@pytest.mark.parametrize(
    "cmap, cmap_var", [(plt.get_cmap("Set1"), 10), (plt.get_cmap("viridis"), 300)]
)
def test_manhattan_accepts_cmap_var_above_colormap_size(cmap, cmap_var):
    fig, ax = plt.subplots()

    qqman.manhattan(make_assoc([1, 2, 3]), ax=ax, cmap=cmap, cmap_var=cmap_var)

    assert len(ax.collections) == 3


@pytest.mark.parametrize("plot", [qqman.manhattan, qqman.qqplot])
def test_out_without_extension_in_dotted_directory_saves_png(plot, tmp_path):
    out = tmp_path / "plots.v2" / "figure"
    out.parent.mkdir()

    plot(make_assoc([1, 2]), out=str(out))

    assert out.read_bytes().startswith(b"\x89PNG")


def test_qqplot_without_valid_p_values_raises_clear_error():
    fig, ax = plt.subplots()

    with pytest.raises(Exception, match="No p-values between 0 and 1"):
        qqman.qqplot([0.0, 1.0, np.nan], ax=ax)


@pytest.mark.parametrize("kwargs", [{"color": "blue"}, {"c": "blue"}])
def test_qqplot_point_color_can_be_overridden(kwargs):
    fig, ax = plt.subplots()

    qqman.qqplot(make_assoc([1]), ax=ax, **kwargs)

    np.testing.assert_array_equal(
        ax.collections[0].get_facecolor(), [matplotlib.colors.to_rgba("blue")]
    )


def test_manhattan_warns_about_p_values_it_cannot_plot(capsys):
    df = make_assoc([1, 2])
    df.loc[[3, 7], "P"] = 0.0
    fig, ax = plt.subplots()

    qqman.manhattan(df, ax=ax)

    assert "2 SNPs with P <= 0" in capsys.readouterr().out


@pytest.mark.parametrize(
    "chroms, labels",
    [
        ([3, 1, 2], ["1", "2", "3"]),
        (["2", "X", "10", "1", "MT", "23"], ["1", "2", "10", "23", "MT", "X"]),
        (["chr2", "chrX", "chr10", "chr1"], ["chr1", "chr2", "chr10", "chrX"]),
    ],
)
def test_manhattan_orders_chromosomes_numerically(chroms, labels, tmp_path):
    path = tmp_path / "assoc.txt"
    make_assoc(chroms).to_csv(path, sep=" ", index=False)
    fig, ax = plt.subplots()

    qqman.manhattan(str(path), ax=ax)

    assert [t.get_text() for t in ax.get_xticklabels()] == labels


def test_manhattan_keeps_a_categorical_chromosome_order():
    df = make_assoc(["1", "2", "X"])
    df["CHR"] = pd.Categorical(df["CHR"], categories=["X", "2", "1"])
    fig, ax = plt.subplots()

    qqman.manhattan(df, ax=ax)

    assert [t.get_text() for t in ax.get_xticklabels()] == ["X", "2", "1"]


def test_manhattan_gap_applies_to_small_chromosomes():
    fig, ax = plt.subplots()

    qqman.manhattan(make_assoc([1, 2], per_chrom=50), ax=ax, gap=10)

    chr1, chr2 = (c.get_offsets()[:, 0] for c in ax.collections)
    assert chr2[0] - chr1[-1] == 11
