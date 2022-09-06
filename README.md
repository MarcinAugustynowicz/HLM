## Description
This repository contains scripts used to verify Proposition 7.2 from the paper "Homological Lagrangian monodromy groups for some monotone tori" (https://arxiv.org/pdf/2201.10507.pdf).

The dim.gap script extracts all finite subgroups of GL(n, Z), for n = 4,5,6, picks those that are not mentioned in the proposition, and saves its elements and generators.

Then the python script algorithm.py applies Lemma 4.3 with k=1,2, which results in all groups excluded, hence proving the Proposition 7.2. 

# Requirements

You need to download the RatProbAlgTori package from https://www.math.kyoto-u.ac.jp/~yamasaki/Algorithm/RatProbAlgTori/ and put it in the same folder as the "GAP" folder.

# Usage

~~~
gap GAP/dim.gap
python python/algorithm.py
~~~
