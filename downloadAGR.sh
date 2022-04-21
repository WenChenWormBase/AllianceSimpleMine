cd /home/wen/agrSimpleMine/src/
wget http://download.alliancegenome.org/3.2.0/DISEASE-ALLIANCE/COMBINED/DISEASE-ALLIANCE_COMBINED_37.tsv
wget http://download.alliancegenome.org/3.2.0/EXPRESSION-ALLIANCE/COMBINED/EXPRESSION-ALLIANCE_COMBINED_37.tsv
wget http://download.alliancegenome.org/3.2.0/GENE-DESCRIPTION-TXT/WB/GENE-DESCRIPTION-TXT_WB_9.txt
wget http://download.alliancegenome.org/3.2.0/GENE-DESCRIPTION-TXT/ZFIN/GENE-DESCRIPTION-TXT_ZFIN_9.txt
wget http://download.alliancegenome.org/3.2.0/GENE-DESCRIPTION-TXT/FB/GENE-DESCRIPTION-TXT_FB_9.txt
wget http://download.alliancegenome.org/3.2.0/GENE-DESCRIPTION-TXT/HUMAN/GENE-DESCRIPTION-TXT_HUMAN_9.txt
wget http://download.alliancegenome.org/3.2.0/GENE-DESCRIPTION-TXT/MGI/GENE-DESCRIPTION-TXT_MGI_9.txt
wget http://download.alliancegenome.org/3.2.0/GENE-DESCRIPTION-TXT/RGD/GENE-DESCRIPTION-TXT_RGD_9.txt
wget http://download.alliancegenome.org/3.2.0/GENE-DESCRIPTION-TXT/SGD/GENE-DESCRIPTION-TXT_SGD_9.txt
wget http://download.alliancegenome.org/3.2.0/INTERACTION-MOL/COMBINED/INTERACTION-MOL_COMBINED_12.tar.gz
wget http://download.alliancegenome.org/3.2.0/INTERACTION-GEN/COMBINED/INTERACTION-GEN_COMBINED_16.tar.gz
wget http://download.alliancegenome.org/3.2.0/ORTHOLOGY-ALLIANCE/COMBINED/ORTHOLOGY-ALLIANCE_COMBINED_37.tsv
wget http://download.alliancegenome.org/3.2.0/VCF/WBcel235/VCF_WBcel235_38.vcf
wget http://download.alliancegenome.org/3.2.0/VCF/GRCz11/VCF_GRCz11_38.vcf
wget http://download.alliancegenome.org/3.2.0/VCF/R6/VCF_R6_38.vcf
wget http://download.alliancegenome.org/3.2.0/VCF/GRCm38/VCF_GRCm38_38.vcf
wget http://download.alliancegenome.org/3.2.0/VCF/Rnor60/VCF_Rnor60_41.vcf
wget http://download.alliancegenome.org/3.0.0/GENECROSSREFERENCE/COMBINED/GENECROSSREFERENCE_COMBINED_115.tsv
wget https://www.dropbox.com/s/mpeqhkay0n1arni/Alliance_gene_synonyms_stage_Oct_18_2020.csv?dl=0

gunzip INTERACTION-GEN_COMBINED_16.tar.gz
tar -xvf INTERACTION-GEN_COMBINED_16.tar
gunzip INTERACTION-MOL_COMBINED_12.tar.gz
tar -xvf INTERACTION-MOL_COMBINED_12.tar
cat GENE-DES* > Alliance_gene_description.txt
cat VCF* > Alliance_variants.txt

rm diseaseFile
ln -s /home/wen/agrSimpleMine/src/DISEASE-ALLIANCE_COMBINED_37.tsv   diseaseFile
rm exprFile
ln -s /home/wen/agrSimpleMine/src/EXPRESSION-ALLIANCE_COMBINED_37.tsv exprFile
rm orthoFile
ln -s /home/wen/agrSimpleMine/src/ORTHOLOGY-ALLIANCE_COMBINED_37.tsv orthoFile

