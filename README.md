This is the primary source code repository for the paper [Numerical equality tests for rational maps and signatures of curves](https://arxiv.org/pdf/2005.04783.pdf), which appeared in the proceedings of the ACM Conference [ISSAC 2020](https://dl.acm.org/conference/issac/proceedings).

This implementation is for the [Macaulay2](http://www2.macaulay2.com/Macaulay2/) computer algebra system and is written in the top-language Macaulay2 language. We recommend the latest version 1.16, though 1.15 may also be fine. 

Commented examples illustrating the basic functionality maybe found in "teaser.m2", "fermat.m2", "quadrifolium.m2". These can be run line by line or as scripts.

There are three main routines whose usage is illustrated in these example files---[witnessCollect](https://github.com/timduff35/NumericalSignatures/blob/master/main.m2#L17), equalityTest, and runMonodromy.

Data for the figures in section 4 may be reproduced using the scipts in the experiment_scripts subdirectory. For each m2 file (eg "joint-3pan.m2"), simply start a m2 session in this directory and excecute the lines below the "end" statement. Each experiment takes 30-120 minutes. Each experiment is reproducible, so you should obtain the same data on false negatives that we report. Note however that timing data may differ. Summaries of this data and the last figure may be obtained via the provided R scripts.
