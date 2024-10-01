#!/usr/bin/env python

import sys
import polars as pl
import altair as alt

###########################

def first_line(file):
    with open(file) as f:
        return f.readline().rstrip()

def build_histogram(data: pl.DataFrame, filename: str, max: int):
    histogram = data["frequency"].hist(bins=[x for x in range(0, max + 1)], include_category=False)
    histogram = histogram.select("breakpoint", "count").filter(histogram["breakpoint"] <= max)
    plot = alt.Chart(histogram).mark_bar().encode(
        x=alt.X("breakpoint"),
        y=alt.Y("count").scale(type="symlog"),
    )
    plot.save(filename)


files = sys.argv[1:]

schema = {"id_high_degree": pl.Int64, "id_low_degree": pl.Int64, "k": pl.Int64, "frequency": pl.Int64}
accumulated = pl.DataFrame({}, schema=schema)
summed_max = 0
for f in files:
    data = pl.read_csv(f, comment_prefix="#", schema={"id_high_degree": pl.Int64, "id_low_degree": pl.Int64, "frequency": pl.Int64})
    curr_max = int(first_line(f).split("max=")[1])
    summed_max += curr_max
    k = None
    for i in range(2, 128):
        if f"k{i}" in f:
            k = i
    data_with_k = data.with_columns(k=pl.repeat(k, data.select(pl.len()).item(), dtype=pl.Int64))
    accumulated = pl.concat([accumulated, data_with_k], how="diagonal")
    # build_histogram(data, f + ".histogram.pdf", curr_max)

summed = accumulated.group_by("id_high_degree", "id_low_degree").agg(pl.col("frequency").sum())
# build_histogram(summed, "summed.histogram.pdf", summed_max)
print("id_high_degree,id_low_degree,frequency")
print(f"# max={summed_max}")
for row in summed.iter_rows():
    print(f"{row[0]},{row[1]},{row[2]}")
