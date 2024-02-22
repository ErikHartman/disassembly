# Protein degradation simulation and complexity estimation

- `src`for code.
- `paper`for paper and related files.
- `data`for data.

- `timeline.pptx` is for notes.

## Background

The peptidome is generated as proteases degrade proteins. While proteins generally are degraded to re-use amino acids and regain energy, many peptides also have important biological functions. Importantly, they play a part of the innate immune response by providing antimicrobial avtivity and being immunomodulatory. This makes us think that protein degradation is a highly regulated process, where evolution has resulted in the selection of certain degradation products and paths. We (and others) have shown that degradation patterns differ under certain physiological conditions, such as infection, making us hypothesize that protein degradation is dynamic and is an effect of the environment.

Similarly to how antibody repertoires hone in on certain antibodies to fight an infection, we wonder if the peptidome is regulated to focus on certain functions during the infection. Do we have a general turnover of the proteome in a normal state, which then shifts to highly selective degradation in certain conditions?

In oct 2023, a paper was published to show that assembly theory can be used to quantify selection. We wondered if a similar theory could be created to quantify selection in degradation of the proteome.