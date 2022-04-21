
Alliance SimpleMine Source Data Building Pipeline

1. prepSimpleMineFileStructure.sh
This shell script creates the source data file structure

2. downloadAGR.sh
This shell script download files from the AGR Download page, store them in the src/, and create symbolic links for them.

3. prepInteraction.pl
This perl script work on interaction files downloaded from AGR Download and create the InteraInteraction.csv file for the next step.

4. makeSimpleMineTable.pl
This perl script create the source data for Alliance SimpleMine and populate the whole file structure created by step 1.

5. category_headers
This is the category header for the SimpleMine web script to display search fields.





-----------------------------
Here are the symbolic links for the downloaded AGR data:
lrwxrwxrwx 1 wen wen 57 Oct 20 18:11 desFile -> /home/wen/agrSimpleMine/src/Alliance_gene_description.txt
lrwxrwxrwx 1 wen wen 60 Oct 21 13:47 diseaseFile -> /home/wen/agrSimpleMine/src/DISEASE-ALLIANCE_COMBINED_37.tsv
lrwxrwxrwx 1 wen wen 63 Oct 21 13:47 exprFile -> /home/wen/agrSimpleMine/src/EXPRESSION-ALLIANCE_COMBINED_37.tsv
lrwxrwxrwx 1 wen wen 61 Oct 20 18:10 geneticIntFile -> /home/wen/agrSimpleMine/src/alliance_genetic_interactions.tsv
lrwxrwxrwx 1 wen wen 63 Oct 20 18:11 molIntFile -> /home/wen/agrSimpleMine/src/alliance_molecular_interactions.tsv
lrwxrwxrwx 1 wen wen 62 Oct 21 13:47 orthoFile -> /home/wen/agrSimpleMine/src/ORTHOLOGY-ALLIANCE_COMBINED_37.tsv
lrwxrwxrwx 1 wen wen 77 Oct 21 11:02 synonymFile -> /home/wen/agrSimpleMine/src/Alliance_gene_synonyms_stage_Oct_18_2020.csv?dl=0
lrwxrwxrwx 1 wen wen 49 Oct 20 18:15 varFile -> /home/wen/agrSimpleMine/src/Alliance_variants.txt
lrwxrwxrwx 1 wen wen 63 Oct 20 18:13 xrefFile -> /home/wen/agrSimpleMine/src/GENECROSSREFERENCE_COMBINED_115.tsv

