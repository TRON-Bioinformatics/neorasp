import pandas as pd
import glob
import os
import sys
import plotly.express as px
from plotly.offline import plot

search_path = sys.argv[1]
# Find all benchmark files
bench_files = glob.glob(os.path.join(search_path,"results/*/benchmark/*_bench.txt"))

records = []
for f in bench_files:
    rule = os.path.basename(f).replace("_bench.txt", "")

    # Remove scatter indices like _00, _01 etc.
    base_rule = rule.rsplit("_", 1)[0] if rule.split("_")[-1].isdigit() else rule

    df_file = pd.read_csv(f, sep="\t")
    df_file["rule"] = base_rule
    df_file["file"] = f
    records.append(df_file)

df = pd.concat(records, ignore_index=True)

### Boxplots ###
fig_time = px.box(
    df, x="rule", y="s", points="all",
    title="Job duration",
    labels={"s": "Execution time (s)", "rule": "Rule"}
)

fig_rss = px.box(
    df, x="rule", y="max_rss", points="all",
    title="Maximal ram usage per rule",
    labels={"max_rss": "Max RSS (kB)", "rule": "Rule"}
)

fig_cpu = px.box(
    df, x="rule", y="cpu_time", points="all",
    title="CPU time",
    labels={"cpu_time": "CPU time (s)", "rule": "Rule"}
)

fig_io_write = px.box(
    df, x="rule", y="io_out", points="all",
    title="I/O write",
    labels={"io_out": "Written bytes", "rule": "Rule"}
)

fig_io_read = px.box(
    df, x="rule", y="io_in", points="all",
    title="I/O reading",
    labels={"io_in": "Read bytes", "rule": "Rule"}
)



html_content = f"""
<h1>NeoRasp  Report</h1>
<p>Resource usage of rules (Boxplots)</p>
<h2>Execution and Job Duration</h2>
{plot(fig_time, output_type='div')}
<h2>RAM</h2>
{plot(fig_rss, output_type='div')}
<h2>CPU</h2>
{plot(fig_cpu, output_type='div')}
<h2>I/O usage</h2>
{plot(fig_io_read, output_type='div')}
{plot(fig_io_write, output_type='div')}
"""

outfile = "snakemake_benchmark_boxplot_report.html"
with open(outfile, "w") as rep:
    rep.write(html_content)

