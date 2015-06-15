# The HistoneDB 2.0: a phylogeny-based resource for all histone proteins and their variants #

Eli J. Draizen<sup>1</sup>, Alexey K. Shaytan<sup>1</sup>, Leonardo Marino-Ramirez<sup>1</sup>, Paul B. Talbert<sup>2</sup>, David Landsman<sup>1</sup>, Anna R. Panchenko<sup>1</sup>

<sup>1</sup>Computational Biology Branch, National Center for Biotechnology Information, National Library of Medicine, National Institutes of Health, 8600 Rockville Pike, MSC 6075, Bethesda, MD 20894-6075, USA,
<sup>2</sup>Howard Hughes Medical Institute, Basic Sciences Division, Fred Hutchinson Cancer Research Center, Seattle, WA 98109, USA.

#### Contact ####
Eli Draizen (eli.draizen@nih.gov) or Anna Panchenkeno (panch@ncbi.nlm.nih.gov)

### Abstract ###
Since the initial creation of the Histone Database, several variants, or non-allelic protein isoforms, for each core histone have gained new interest. These variants have have a vast number of functions and play important roles in nucleosome dynamics. We need the Histone Database to account for these new variants. Here, The Histone Database has been renamed HistoneDB and rebuilt to be organized and searchable by variants and taxonomy. We created profile-HMMs for each variant to search the nr database. Sequences were added into the HistoneDB if they were found to have 95% specificity. HistoneDB 2.0 is built with the Django Python Web Framework and the source is available at GitHub: http://github.com/edraizen/HistoneDB.

Database URL: http://www.ncbi.nlm.nih.gov/projects/HistoneDB 

## Introduction ##
The traditional role of histone proteins is in their formation of the the nucleosome. Histones were previosly thought to only have a sturcutral role, compacting and organizing the DNA to fit inside the nucleus. However, as we uncover more information about histones, we have learned they are necessary to controlling and regulating gene expression, cell reprogramming, and much more.

Fig 1. Chromatin and nucleosome structure. Modified from From Lewins Genes X & Molecular Cell Biology, Lodish H et al. 
 
The nucleosome core particle is formed from an octamer of four core histone types: H2A, H2B, H3, and H4 that all share the same histone fold, but contain less than 25% sequence similarity. The nucleosome is formed by two copies of each histone, with H2A-H2B forming a dimer (2x) and H3-H4 form a tetramer. This octamer of histone proteins wraps around 189 basepairs of DNA 1.67 times. There is also a linker histone H1, which sits above the linker DNA of the nucleosome, but does not share the histone fold. \ref{Shaytan2015}

Fig 1. The Nucleosome core particle. H2A is yellow, H2B is red, H3 is blue, and H4 is green. H1 is not shown. Figure created by \ref{Shaytan2015} 

A histone variant is defined as a non-allelic protein isoforms that forms monophyletic clade. Each varaint  has a specific structure and function that is different than the canonical. It is important to note that histone variants can become post translationally modifications (PTMs), but this database does not include PTMs at this time. 

Another important aspect of this database is to promote the new nomenclature for histone variants defined by Talbert et al based on phylogeny. If a clade is further separated into smaller clades, those are defined as genes and splice isoforms. Each branch point becomes a period in the new notation.

#### H2A ####
There are six structurally distinct variants for H2A. H2A.X is the most common, with notable sequence motif ‘SQ(E/I)Y’. H2A.X is known to be involved with DNA damage response, chromatin remodelling, and X chromosome inactivation in somatic cells. H2A.X, however, is the only variant that is not monophyletic; its has evolved several times, but each version has a similar structure and function. Currently only the final motif SQEY has a PDB.

H2A.Z is known for transcription regulation, DNA repair, suppression of antisense RNA, and Polymerase II recruitment. Notabel features for H2A.Z include a large hydrophobic patch, a sequence motif of ‘DEELD,’ a one amino acid insertion in loop1,and a one amino acid deletion in the docking.

macroH2A contains a histone histone fold, c-term macro domain, can bind ADP, lost in several lineages including Caenorhabditis and Drosophila \ref{Talbert2012}, X-inactivation, Positive or negative transcriptional regulation. PDBs of each domain, but linker is too flexible to be crystalized

H2A.B is a rapidly evolving **B**arr body deficient varaint found in mammals. It is known for its involvment with spermiogensis. H2A.B has large expansions and a shortened docking domain, wrapping less DNA. It is colsely related to H2A.L and H2A.M.

H2A.L: mammals, mammalian spermiogenesis, rapidly evolving, large expansions, shortened docking domain so wraps less DNA, forms subnucleosomal particles with TS H2B.1. Homology model available?

H2A.M is the newest mamlain-speific varaint to date. It binds to huntingtin protein M.

other less studied H2A variants include H2A.J, which is very similar to canonical, TS H2A.1, which is testis-specific, and H2A.Q, which had been identified in Oikopluera. There are also species-specific variants, H2A.1 through H2A.10, however the number does not imply homology

#### H2B ####
TS H2B.1, PDB available
H2B.W: Spermiogenesis, Telomere associated functions in sperm, found in Spermatogenic cells
H2B.Z
subH2B: Spermiogenesis, found in Nucleus / subacromosome of spermatozoa, Has a bipartite Nuclear localization signal:

#### H3 ####
cenH3, or centromeric H3, is found when the nucleosome is forming the centromere. cenH3 is known to bind the protein CENP-C. Currently, the cenH3s do not form a monophyletic clade, but many beileve this is due to limited phylogentic techniques, and do actually form a monophyltic clade. While cenH3s are more of a functional class, they do contain strucutral features such as an extended loop1.

H3.3 is another well studied variant that has diverged multiple times. The origin of canonical H3s appears to me to have occurred separately in plants, animals, and other lineages. I would not say they are species-specific in most cases; nearly all animals have the same H3.2 and H3.3 proteins. H3.2s and H3.3s are identical between flies, humans and most other animals

H3.3 PDN Available
H3.Y: 
H3.1:
H3.2:

Canonical H3s and H3.3s have diverged multiple times. " -Talbert

#### H4 ####
mostly single form of H4, there are some species specifc varaints. PDB
H4.1
H4.2
H4.V: found in Kinetoplastids

#### H1 ####
H1

Each variant has distinguishable features that can be captured by a Hidden Markov Models. All of the sequences from the non-redundant (nr) database have been classified by our variant models and added into the HistoneDB so it can be easily searchable by variants.

 

However, if you are still uncomfortable with new naming scheme, you can search the database using the old naming schemes. We will tell you the new correct name, so you learn to adopt the new standard.

While this update does not deal with histone-like proteins in archaea and bacteria, you can still access the sequences collected for these types developed for the previous version of the database.

## Database and Software ##
### Data sources ###
We collected sequences with a broad phylogenetic representation for each variant of H2A, H2B, H3, H4, and H1, which were manually curated using blast \ref{Talbert2012}. Starting with known variants, e.g. from Human or Arabidopsis, we used tblastn and blastp with psi-blast to search nr for similar variants. Model organisms were favored as well as non-partial sequences, and sequences with the correct splice isoform. In most cases, the taxonomic limiting box was used to target certain groups or exclude the abundance of identical vertebrate proteins. For more divergent eukaryotic groups like Alveolates and Excavates, an identified histone within the group was used to search for others. This strategy works well with closely related organisms, but in other cases, it simply increases the divergence between proteins that have diverged from a "consensus" sequence in different amino acids. For example, when looking for cenH3s, we used a canonicalH3, rather than a cenH3, to search against a targeted taxonomic group and indentified hits that had about 50% identity. We also made sure that the hits contained notable sequence features that are known for each variant. These sequences were also used to make phylogentic trees in \ref{Talbert2012}. 

Previous variant names were manually extracted from \ref{Talbert2012}.

We also used the core and linker histone multiple sequence from the previous Histone Database to find canonical histone variants. 

### Histone variant identification ###
Once each variant was in its respective FASTA file, each variant was separately to create seed alignments for each variant. These were checked manually to make sure they had a wide taxonomic distribution and no large insertions in the core histone fold regions. These seed sequences were then used to train profile Hidden Markov Models, using HMMER 3.1b2 hmmbuild, creating one profile-HMM for each variant. Next, all of the variant models were combined into one file and pressed using HMMER 3.1b2 hmmpress. 

Further we defined cutoffs scores for each variant based on 95% specificity. Positive training sets were defined as the seed for the given variant and negative examples we defined as a combination of all seeds except the given varaint \fig{H2AZ_cutoff). Iterating over thresholds from 0 to 1000, we calculated the specificty at each threshold, interpolating the curve to predict all speficifty values. We chose the cutoff score when the specifity was equal to 95%. Please see supplematary information for details about evaluting each varaint model.

Finally, we used the combined HMM to search all of the NCBI non-redundant (nr) database using HMMER 3.1b2. To classify each sequence, we took the maxium score from a varaint model that was also above the threshold for that variant. In other words, we iterated over the hmmsearch output from each variant, saving the sequence into the database and either classifying the sequence if the sequence has not been seen before if the score against the current variant model was greater than the threshold for the varaint or re-classifying the sequence if the sequence was already classified, the score was greater than the previous classified variant, and the score was greater than the threshold for the current variant. To see a workflow, please see \ref{classification}. If the sequence was not able to be classified, the sequence was stored with variant "Unknown".

To find canonical histones we built core histone profile-HMMs using the alignments from the original Histone Database to search nr again. Since these alingments are for each core histone, the scores against this model are much larger than they would be against each variant model. We avoided this problem by classying all the variants first, and then canoncal after, using the same testting and classification method, but ignoring sequences that were already classified.

### Software ###
The HistoneDB 2.0 is now written in Django, a high-level Python Web Framework, with a MySQL backend. The project has two applications, ‘browse,’ and ‘djangophylcore.’ Browse contains the HistoneDB models (equivalent to database tables), views (python functions to display each page), and templates (HTML files). Djangophylocore is a previously developed django application developed at University Montpellier to store the NCBI taxonomy database in a Django relational database using an algorithm similar to Modified Preorder Tree Traversal. The scheme for the HistoneDB Django database can be seen in \ref{django_schema}.

Figure x. The HistoneDB database schema. 

The layout is based on Twitter Bootstrap, with four important pages: Main browse of all core histone types, Core histone type browse of a single core histone, Variant browse, and All Sequences/Search.

#### Main browse ####
The frontpage is a general browse page, where you can choose one of the four core histone types, by selecting a color coded 3D model for each variant created from PDB 1AOI (Shaytan, 2015). You can also select archaeal and bacterial sequences and sequences that contains structures in the PDB, which were all curated for the last version of the database. 

#### Core Histone browse ####
After you choose a core histone type, you are redirected to a page for the individual histone type with four tabs. The first tab contains the phylogenetic tree organized by variants using jsPhyloSVG. These trees have been created by aligning all of the seed sequences for a given histone type using MUSCLE v3.8.31 and neighbor-joining in CLUSTALW 2.1 to create the trees. Th trees were then converted to PhyloXML using BioPython and were edited to add features in jsPhyloSVG. Selecting a variant goes to a new variant page.

The next tab contains all organisms that contain the core histone type as a zoomable D3 sunburst. These have been created using the Djangophylore parent attribute for each sequence following back the root, only allowing taxonomies with a rank and stopping at the rank 'order.' Selecting a node will zoom into a new sunburst with the selected node in the center and it’s children as leaves. When you have finished drilling down the sunburst, you hit ‘View Sequences,’ which will take you to the ‘All Sequences’ page that will display sequences that have the core type you originally selected and are within the selected phyla.

The third tab contains all sequences of the selected core histone type. This is based off the ‘All Sequences’ page.

The fourth tab contains a general seed alignment for the core histone fold. These were previously calculated for the initial Histone Database.

#### Variant Browse ####
The variant browse page starts off with the name of the variant, a description of its proposed functions, and previous names the variant has been referred to. Next you have three tabs, with the first being an organism tab that displays a D3 sunburst like the core histone browse page, but this one only contains organism that have this variant. Again, after you have zoomed into the desired organism, you can select ‘View Sequences,’ which will take you to the ‘All Sequences’ page that will display sequences that have the core type you originally selected, the selected variant, and are within the selected phyla.

The next tab contains all of the sequences for that variant, with the same functionality as the core histone browse page.

The final tab consists of the multiple sequence alignment for the specific variant, which have been created using the sequences identified by Talbert. This is also made with the BioJS msa module.

#### All Sequences ####
You are presented with the GI, Variant, Gene, Splice, Species, Header, and Score for each sequence using wenzhixin’s bootstrap-table. You can then select one or more sequences to see its multiple sequence alignments with secondary structure annotations, scores for the sequence in other models, view in Entrez, and download all. The Multiple Sequence Alignment is made from the BioJS msa package. Viewing scores of the sequence from other models is important to see if there was a classification error, or there may be other significant variants but we don’t know for sure, e.g. the sequence with GI xxxx from Drosophila is classified at H2A.Z, but it contains the motif “SQEY,” which is common in H2A.X.

#### Search ####
There are two types of search, ‘filter’ and ‘simple search’. Filter is the standard search, where you presented by database column names to filter by. You can also use the filter modal to search for high scoring sequences in other variants. This is used in Advanced Search and the All Sequences Advanced Filter, where the Adanced Search takes you to a new seach page with results and the Advanced Filter filters the current view, which is saved in a sission variable. Simple search is when you have a single textbox that filters the entire database, which can be seen on the navigation bar and the All Sequences search textbox. The simple search searches in order GI number, core histone type, variant, old name schemes for variant, taxonomy, sequence headers, and finally sequence motifs, stopping if a match is made. If the search was from the navigation bar and the result was a core histone type or variant, the page is redirected to its respective browse page. 

## Results and Discussion ##
We were able to classify most histone sequences by variant in the nr database.For the varaints we were able to classify, we were able to show.

H2A | | H2B | | H3 | | H4 | | H1
----|-|-----|-|----|-|----|-|----
H2A.B | 121 | H2B.1 | 7412 | cenH3 | 256 | canonicalH4 | 2990|  
H2A.L | 165 | H2B.W | 216 | H3.1 | 2588 | | |
H2A.M | 72 | H2B.Z | 123 | H3.2 | 1847 | | |
H2A.X | 7690 | subH2B | 86 | H3.3 | 14 | | |
H2A.Z | 2324 | canonicalH2B | 3301 | H3.Y | 21 |
macroH2A | 1135 | | | canonicalH3 | 36717 | | | |
canonicalH2A | 0 | | | | |

Compared to the previous Histone Database, we have found **X** more variants, and **Y** more sequences.  


Our models, however, are not able to distinguish between paralogous genes and splice isoforms, which we we call gene and splice isofform in the database.

## Conclusion ##
The HistoneDB has been redesigned to organize the histones in nr by variant; provide reference alignments for each variant; understand how the variants evolved; find orthologs of variants in other species; and finally, promote the new histone varaint nomenclature \ref{Talbert2012}. 

This database will allow chromatin researchers to find sequences for each variant from different organsms which can be used for 3D modelling. A rising field is the study of nucleosome dynamics through molecuar dynamics simulation. Since now we have an annoted set of histone varaints, researchers can extend this study to see what effects of different varaints have nucleosme dynamics. It would also be intrguing to understand which varaints dimerize together to form the complete octomeric structre of the nucleosome. 

The database can be easily extended to include new variants, and rebuilt if the nr database changes.

## Funding ##
This work was supported by the Intramural Research Program of the National Library of Medicine, NIH. ED is supported by the Oak Ridge Institute for Science and Education. AS is supported by the US-Russia Collaboration in the Biomedical Sciences NIH visiting fellows program.

## Acknowledgements ##
We would like to thank Franco Simonetti for discussions about Django and the IT support at NCBI for troublshotting server ralted issues.

ED wrote the manuscript. PT collected the sequences and wrote the 'data sources' section via summarized email communication. AS, LM, DL, and AP supervised the project. 

Conflict of interest~. None declared.

## References ##


