---
title: Disassembly theory and its appliction to protein degradation

output: word_document
bibliography: bibliography.bib
---

Assembly theory (AT) has recently been introduced as a novel field which aims to quantify selection and evolution. In short, it is a theoretical framework which looks at two properties of an object: its copy number and its assembly index (an assembly index is the number of steps on a minimal path producing the object). 

![](/paper/img/at3.png)
*From Fig. 1 in Cronin et al. An example of assembly through recursive steps.*

The [paper by Cronin](https://www.nature.com/articles/s41586-023-06600-9) [@assembly] introducing AT was published in Nature in Nov. 2023 and got a lot of attention (~200k reads in a month).

Here are some important quotes from that paper:

- *"The concept of an object in AT is simple and rigorously defined. An object is finite, is distinguishable, persists over time and is breakable such that the set of constraints to construct it from elementary building blocks is quantifiable."*

- *"In AT, we recognize that the smallest unit of matter is typically defined by the limits of observational measurements and may not itself be fundamental. A more universal concept is to treat objects as anything that can be broken and built. This allows us to naturally account for the emergent objects produced by evolution and selection as fundamental to the theory."*

- *"The concept of copy number is of foundational importance in defining a theory that accounts for selection. The more complex a given object, the less likely an identical copy can exist without selection of some information-driven mechanism that generates that object."*

![](/paper/img/at2.png)
*From Fig. 2 in Cronin et al. A high copy number in complex molecules suggests selection.*

![](/paper/img/at1.png)
*From Fig. 5 in Cronin et al. Assembly without and with selection.*


There is a lot of potential utility in AT, the most apparent is the ability to quantify selection through the assembly equation:

$A = \sum_{i=1}^Ne^{a_i}(\frac{n_i-1}{N_T})$

where $A$ is the assembly, $a_i$ is the assembly index, $n_i$, the copy number of object _i_ and $N_T$ the total number of objects in the ensemble.

This could be used to quantify selection on e.g., a different planet, by analyzing the molecular space with a mass spectrometer.

## Disassembly
I propose a subfield of AT which deals with the disassembly of fundamental objects, e.g., the disassembly of a protein into peptides, named disassembly theory (DT).

DT has many of the same aspects of AT, except that selection occurs through the selective disassembly of larger objects, essentially doing AT backwards. If we have a large quantity of similar degradation products with high copy number, this implies that there is some selection for the degradation product. This could be used to e.g. quantify the selection processes which has given rise to a peptidome.

The fundamental object is in this case something which cannot be combined to make new things, or alternatively, things which have no parent-thing (i.e. is not a part of a greater thing). This is e.g. a protein, since two proteins cannot be combined to make a new protein, and no protein sequence is a sub-sequence of a larger protein.

The main difference between AT and DT is the assembly index/disassembly index.

To make a complex object you need many steps, but to generate a small object from a larger one you may just need 1 or 2 steps. E.g., a peptide has two cut-sites (if non-terminal), and any peptide can theoretically be generated from a protein using just two steps.*

_*if we don't consider the 3D structure and cut site availability of the protein._

Degradation may also be sequential, where a larger degradation product is subsequentially degraded into a smaller object. This resembles the notion of assembly index in AT. However, this cannot be proven without evidence for parent objects.

To encompass these phenomena, we can modify the assembly equation into the disassembly equation. First we define the disassembly index, $d$, as the length of the number of possible paths from one object to the fundamental object in the ensemble.

Since we don't know what path has been taken, we can use either:

- the longest path available in the data from one object to the parent object
$\delta_i = max(d_i)$
- or the mean path from the object to the parent object.
$\delta_i = \frac{1}{N} \sum_{i=1}^{N}d_i$


The disassembly, $D$, can then be calculated as follows:

$D = \sum_{i=1}^N e^{\delta_i}(\frac{n_i-1}{N_T})$

![](/paper/img/sketch.png) 

I hypothesize that if we have selection in the peptidome, we will see evidence of this in the disassembly index, and this could be used to quantify that selection. I think that e.g. a bacteria with a lot of impact on the proteolytic environment will generate a peptidome with larger $D$ than a bacteria with less impact. Controls can indicate baseline selection by endogenous proteases.

## Experiments

Theoretically, I think disassembly is sound and has application on peptidomics. To test this I will start with some experiments on synthetic data. 

I will then apply it to the real wound data.

This could end up as a short 3 figure paper. 2 theoretical and 1 application. Strike the iron while it's hot. The original paper has 4 citation since release (1 month ago).

