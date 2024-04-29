This directory contains detailed information on environment packages, including versions, for different environments (conda envs: yml files, guix env: list of packages with their versions).

[You can recreate those conda environments from the yml files.](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)

Please consider the following notes.

## scib_v1.0.1_min
Manual fix to line 594 of scib/metrics/lisi.py
Changed:
```print(i[1] + " has not enough neighbors.")```
To:
```print(str(i[1]) + " has not enough neighbors.")```

## dynamic_LIAM_challenge_reproducibility
In addition to the packages listed in the yml file, we installed the package liam_NeurIPS2021_challenge_reproducibility from a local clone of the legacy software [GitHub repository](https://github.com/ohlerlab/liam_challenge_reproducibility), which made changes dynamically available. 
```pip install -e liam_challenge_reproducibility/```

We trained the following models with these respective versions of the liam_NeurIPS2021_challenge_reproducibility package:
- unscaled models + single modality models: v0.0.0
- scaled models: v0.0.1
- mosaic models: v0.0.2

Using [liam v0.1.1](https://github.com/ohlerlab/liam) under active development will give equivalent results as we introduced only downward compatible changes. 