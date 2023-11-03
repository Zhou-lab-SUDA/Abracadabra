# Abracadabra

***A****cinetobacter* ***B****aumannii* **R**esistance **A**nalyse, **C**gMLST **A**ssignment and **D**ecoding **A**ppliaction for **B**ioinformatic **R**esearch **A**pparatus

**Abracadabra**, a spell of ancient origin, is our key to unlocking the enigmatic power concealed within the genomes of *Acinetobacter baumannii*.

In this arcane rite, users must offer their solemn genome sacrifices, whispering their innermost desires. From **species divination** to the intricate arts of **dcgMLST/SNP decoding analysis** and the quest for **resistance genes**, Abracadabra shall bestow upon you the mystical blessings of science.

We envision a future where Abracadabra, like a magical artifact, shall unveil new prophetic abilities to **foresee plasmids**, unraveling further scientific wonders.

# Installation
Abracadabra spell was devoloped and tested in Python 3.9.0, and requires a several modules:
~~~~~~~~~~
numba
numpy
pandas
biopython
pyarrow
fastparquet
argparse
fastANI
ncbi-blast+
diamond
~~~~~~~~~~

Enjoy a coffee while installing these packages using commands below:
~~~~~~~~~~
pip install click numba numpy pandas biopython pyarrow fastparquet fastani
sudo apt install -y ncbi-blast+ diamond
or
conda install -c conda-forge biopython numba numpy pandas pyarrow fastparquet
conda install -c bio-conda blast diamond
FYI: some version of pandas might not be support DataFrame.append() any more.
~~~~~~~~~~

# Quick Start (with examples
##
