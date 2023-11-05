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
samtools
bcftools
minimap2
~~~~~~~~~~

Enjoy a cup of coffee while installing these packages using commands below:
~~~~~~~~~~
pip install numba numpy pandas biopython pyarrow fastparquet fastani minimap2
sudo apt install -y ncbi-blast+ diamond
or
conda install -c conda-forge biopython numba numpy pandas pyarrow fastparquet
conda install -c bio-conda blast diamond samtools bcftools minimap2
FYI: some version of pandas might not be support DataFrame.append() any more.
~~~~~~~~~~

# Quick Start (with examples)

To begin your journey, present your genomic tribute to Abracadabra. You may offer a single genome for analysis using the `-q/--query` parameter:
~~~~~~~~~~
Abracadabra -q [input.genome]
~~~~~~~~~~

Or, if you have a collection of genomes, inscribe their paths in a sacred scroll and offer it to the spell using the `-ql/--query_list` parameter:
~~~~~~~~~~
Abracadabra -ql [genome_list.txt]
~~~~~~~~~~

In the realm of sorcery, Abracadabra empowers conjurers with three mystical abilities:

1. Species Identification (Beware, for those who dare to wield Abracadabra must ensure their intent is sane. This spell is designed exclusively for the genomes of *Acinetobacter baumannii*.)
2. Resistance Analysis (Abracadabra harnesses the wisdom of AMRFinder Plus and VFDB databases to predict resistance genes and virulence factors in the genome. It also unveils the enigmatic realm of resistance mutations that manifest in the *Acinetobacter baumannii* lineage.)
3. DcgMLST Analysis (A new and distributed cgMLST system awaits *Acinetobacter baumannii* genomes, capable of scrutinizing 1030 allelic loci. It discerns differences, confirms the genome's phylogeny, and ventures to unmask the mighty empire of ESL, a superlineage invading the realm of humankind. If the uploaded genome belongs to the lineage of ESL, Abracadabra endeavors to identify its specific lineage or variant. For more information about dcgMLST, follow great magician, Zhemin Zhou, Lord of [Enterobase](enterobase.warwick.ac.uk/))

The path to mastery lies in your ability to choose the right incantation. To wield Abracadabra's powers, you may use the `-t/--tasks` parameter, which accepts the following mystical commands:

- `s` to seek the wisdom of species identification.
- `r` to summon the knowledge of resistance and virulence gene prediction.
- `c` to delve into the arcane realm of dcgMLST analysis.
- you can also combine the commands together, like `src` for the all analysis or `rs` to exclude the dcgMLST.

Beware, for any conjurer who utters a foolish command shall find Abracadabra deaf to their pleas. Choose your path wisely, and the magic shall reveal its secrets.

## First spell for species indentification
~~~~~~~~~~
Abracadabra -q [input.genome] -t s
~~~~~~~~~~
With this mystical spell, conjurer shall unravel the secrets of their genome and its species against the reference genome of **K09-14**(reference genome for *Acinetobacter baumannii*). The incantation operates by measuring the differences between your offering and K09-14, using a threshold of 95% average nucleotide identity. In the realm of species *Acinetobacter baumannii*, values between 97-99% are often the norm, as **K09-14** may not be the regular archetype.

Example output:
~~~~~~~~~~
2023-11-05 14:13:11     ---In the realm of magic, species identification begins---

2023-11-05 14:13:12     The offering genome is revealed as Acinetobacter baumannii

2023-11-05 14:13:12     Average Nucleotide Identity: 97.749

2023-11-05 14:13:12     ---Magic fades away leaving knowledge for Acinetobacter baumannii---
~~~~~~~~~~

## Second spell for ARGs and virulent prediction
~~~~~~~~~~
Abracadabra -q [input.genome] -t r
~~~~~~~~~~
With this mystical spell, conjurer shall embark on a quest to uncover the concealed knowledge of resistance genes, virulence factors, and the enigmatic mutations that give rise to alterations in resistance phenotypes. The spell employs the wisdom of [AMRFinder Plus]() and [VFDB]() databases as its guiding lights. Yet, it does not stop there. The spell extends its gaze to certain  genes and their homologs and certain antibiotics, for their relationship have revealed. For example, the tet(B) gene is believed to  mediate resistance to minocycline, and the ADC-30/73 gene family is believed to be associated with resistance to sulbactam.

Example output:
|contig|accession_number|identity|coverage|mismatch|gap|start|end|evalue|bit_score|gene|function|function_category|description|
|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
|test1|WP_000052512.1|100.0|100.0|0|0|1400750|1402225|0.0|2168|msr(E)|MACROLIDE|MACROLIDE|ABC-F_type_ribosomal_protection_protein_Msr(E)|

## Final spell for dcgMLST/HierCC assignment
~~~~~~~~~~
Abracadabra -q [input.genome] -t c
~~~~~~~~~~
With this powerful spell, conjurer shall embark on a quest to discover your genomic offering's place within the intricate tapestry of *Acinetobacter baumannii*. It is guided by our newly established dcgMLST system. The spell shall unveil the differences between your offering and others, using the threshold of HC1030(allowing no more than 1030 allelic differences between any two genomes in one cluster).

The command shall yield a dcgMLST profile, a parchment bearing the MD5 imprints of core genes, suitable for exchange or submission. The lineage classification shall unfold the secrets of the phylogeny of *Acinetobacter baumannii*. We have constructed the following diagram depicting the phylogenetic of *Acinetobacter baumannii*, you can see some groups considered to be International Clones (IC).
![dcgMLST](https://github.com/Zhou-lab-SUDA/Abracadabra/assets/88537949/4a31a9f4-3322-46ee-99e1-e14bfa12b35e)

For those genomes are deemed as **HC1030_1** or **HC1030_2**, the spell shall further unlock the mysteries of **ESL**(as described below).

Example output:
~~~~~~~~~~
2023-11-05 11:01:22     ---In the realm of magic, cgMLST shows the deep connection inside the bacteria---

2023-11-05 11:02:30     HC1030: 2       Most closed genome: GCF_000805055       Alleles Differed: 0

2023-11-05 11:02:30     ***Behold, the malevolent genome of the ESL has been unveiled!***

                        Let the ancient art of SNP analysis unveil its lineage and variant!

2023-11-05 11:02:41     ESL_lineage: 2.5        ESL_variant: 2.5.2

2023-11-05 11:02:41     ---Writing cgMLST---

2023-11-05 11:02:41     ---Magic fades away leaving knowledge for Acinetobacter baumannii---
~~~~~~~~~~

## ESL definition
HC1030_1 and HC1030_2 formed a monophyletic cluster that consists of a total of 10,820 (69%) strains from 79 countries with highest AMR genes. We therefore denoted the clade as the epidemic super-lineage (ESL) thereafter and focused on its population dynamics. 

## ESL SNP barcode

For those whose offerings are embraced by the lineage of HC1030_1 or HC1030_2, the realm of ESL shall open its doors. We have crafted an enchanted genealogy for ESL, a lineage that holds secrets untold. To retrace the lineage and witness its phylogeny, one may employ the magic of [GrapeTree](https://achtman-lab.github.io/GrapeTree/MSTree_holder.html?tree=github.com/Naclist/Li-et-al.-A.-baumanii-data-repo/blob/main/HC1200_1.r4.json).

