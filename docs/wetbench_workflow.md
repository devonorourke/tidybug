# Laboratory work
## 1. Creating the mock community samples
[Jusino et al., 2019](https://doi.org/10.1111/1755-0998.12951) provide extensive documentation regarding how mock community samples were generated. Briefly, DNA extracted from voucher arthropod specimen was amplified using their novel ANML primers. PCR products were cloned into plasmid vectors and Sanger sequenced to match specimen taxonomy with sequence identity. Taxonomic identities were initially assigned by a trained entomologist’s visual identification. More exclusive taxonomic levels that could not be distinguished by the taxonomists were added by manually aligning full length COI Sanger sequence data to NCBI’s nt database. Sequences were required to have at least 98% identity and 92% coverage to be named to Species level. Any instance in which multiple best hits were available resulted in removing the Species classification (this occurred in just one instance). In addition, we removed a single instance in which a mock sample had a single best hit as "Sp. 1ES" despite being unambiguous because other references in our database and on NCBI suggested there are other "Sp." names that did not clearly identify this record as a single species.

## 2. Guano sampling
Guano was sampled passively by both researchers and citizen scientists as part of a broader collection effort spanning April–October of 2015 and 2016 throughout multiples sites in the United States including New Hampshire, Massachusetts, Vermont, Maine, Indiana, Kentucky, and Virginia. Bat guano DNA extraction and sequencing occurred simultaneously with ongoing collections, and we chose to use only a fraction of these samples (~ 30%) for this analysis because the selected batches analyzed were largely from the same areas (1,675 of 1,790) were collected in sites specific to New Hampshire.

Prospective guano sites were identified through prior participation in summer maternity roost counts by the volunteers. Collection materials were provided and shipped to each volunteer as a kit and included nitrile gloves, face masks, alcohol wipes, forceps, 0.31 mil plastic sheets, collection boxes, and microcentrifuge tubes filled with 1 mL storage buffer (3.5M ammonium sulfate, 16.7 mM sodium citrate, 13.3 mM EDTA, pH 5.2). Volunteers were provided with instructions to lay a sheet of plastic beneath the location with the greatest amount of fresh droppings. We also provided on-site consultations with several volunteers and corresponded with participants electronically to ensure that the sampling procedures were consistent across sites. Sampling started within six days, with up to 10 individual pellets being collected each week (i.e., one vial per pellet). A new sheet of plastic was used each week, and volunteers were instructed to clean forceps between each week of sampling while wearing the provided PPE. Samples were shipped back to our lab and stored at -80 °C until DNA extraction.


## 3. Primer design and amplification of COI gene fragment

We used a dual-indexed primer design described in detail at the [Schloss Lab Github repo](https://github.com/SchlossLab/MiSeq_WetLab_SOP/blob/master/MiSeq_WetLab_SOP.md). This design incorporates both the necessary generic elements of Illumina adapters while allowing for customized barcodes and COI-targeted primer sequences. The COI-specific region of primers match that previously described by [Jusino et al., 2019](https://doi.org/10.1111/1755-0998.12951), resulting in the following generic primer pairs where `<i5>` and `<i7>` represent a unique 8 bp barcode:

```
ANMLxf: 5’-AATGATACGGCGACCACCGAGATCTACAC<i5>TATGGTAATTCGGGTCAACAAATCATAAAGATATTGG-3’)  
ANMLxr: 5’-CAAGCAGAAGACGGCATACGAGAT<i7>AGTCAGTCAGCCGGWACTAATCAATTTCCAAATCC-3’)
```

> Lists for all possible primer pairs with barcodes are listed in the [metadata directory of this project's Github repo](https://github.com/devonorourke/tidybug/blob/master/data/metadata/primerpairs.txt).

We used 25 µL reactions with 10 µL of extracted DNA, 1 µL each of 10 mM forward and reverse primer pairs, and 13 µL of AccuStart II PCR SuperMix (Quanta BioSciences, Gaithersburg, MD, USA). Mock community samples were amplified using a distinct primer pair that was not used on any guano samples for any libraries in this study. Furthermore, mock samples were amplified in a separate reaction from guano samples to avoid any potential cross contamination. For the mock samples, 1 µL of template DNA (10 ng/uL) with 9 µL of nuclease free water was used in lieu of 10 µL of guano extract.

Reaction conditions consisted of:
- an initial denaturation step at 95 °C for 2 min,
- followed by 30 cycles of:
  - denaturation at 95 °C for 20 s,
  - annealing at 50 °C for 15 s, and
  - extension at 72 °C for 60 s;
- and a final extension at 72 °C for 10 min.  

PCR products were quantified using a PicoGreen assay (Invitrogen, Carlsbad, CA, USA) and Tecan plate reader with excitation and emission wavelengths of 480 nm and 520 nm, respectively (Tecan Group, Männedorf, Switzerland). Samples were pooled in approximately equimolar ratios; because there were hundreds of samples pooled in a single library, the exact expected equimolar volumes were rounded to the nearest bin (bin sizes set at 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10, 15, or 20 µL). Extraction blanks or guano samples with detectable DNA that exceeded the maximum pool bin volume were included at a fixed volume of 20 µL. The initial pool volume was reduced with a vacuum concentrator to approximately 2 mL and was cleaned with a QIAquick PCR purification kit (Qiagen, Hilden, Germany); libraries were eluted in 30 µL elution buffer. Libraries were quantified with a Qubit High Sensitivity assay (Thermo Fisher Scientific, Waltham, MA, USA) and fragment sizes were analyzed using TapeStation D1000 ScreenTape (Agilent Technologies, Santa Clara, CA, USA). Note that a single mock community sample was included in each of the four libraries sequenced.

Libraries were submitted to Northern Arizona University and sequenced using an Illumina MiSeq platform (Illumina, San Diego, CA, USA) using v3 chemistry with 600 cycles of 2x300 bp paired-end read lengths. Raw sequence reads available at NCBI [BioProject PRJNA518082](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA518082).
