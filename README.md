# ARR2

Supporting code for the paper "The ARR2 prior: flexible predictive prior definition for Bayesian auto-regressions" by David Kohns, Noa Kallioinen, Yann McLatchie and Aki Vehtari. [Preprint available on arXiv](https://arxiv.org/abs/2405.19920).

Experiments are written in R and Stan and uses snakemake to run.

To run the experiments, enter the snakemake directory and run a rule such as:

```bash
snakemake ar_estim_all # AR estimation experiments
snakemake ar_lfo_all # AR predictive performance experiments
snakemake arx_estim_all # ARX estimation experiments
snakemake arx_lfo_all # ARX predictive performance experiments
snakemake ltx_estim_all # LTX estimation experiments
snakemake ltx_lfo_all # LTX predictive performance experiments
```