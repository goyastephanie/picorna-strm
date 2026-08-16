# Making a Kraken2 human database with hg38 genome + 18s rRNA

## 1. prep
1. download and install [Kraken2](https://github.com/DerrickWood/kraken2)
2. download the [18S Ribosomal N1 Gene sequence](https://rnacentral.org/rna/URS0000726FAB/9606) (DNA)  
    1. on the RNACentral page, select `Copy as DNA` near the bottom of the page.
    2. save this text as a `.fa` file
    3. edit the FASTA file with the taxon ID (required for Kraken2 to add it to database). Here's an example using the current *Homo sapiens* TaxonID 9606:  
    `18sR1rRNA|kraken:taxid|9606`  

>[!IMPORTANT]
>Make sure the text within the pipes `|` are exactly the same as above, and that the taxonID is entered as above.

## 2. creating the initial human database

1. In your chosen work folder, run `kraken2-build --download-taxonomy --db $DBNAME` where `$DBNAME` is the name of the database you wish to create (keep it the same for all ensuing commands).
2. Once the above is done, run `kraken2-build --download-library human --db $DBNAME`.
3. Once the above is finished, *double check that your 18s rRNA fasta contains the appropriate header* and run the following command to add it to the database: `kraken2-build --add-to-library $FASTA_FILE --db $DBNAME`
4. To compile your database to make it usable with Kraken2, run `kraken2-build --build --db $DBNAME`. 

>[!NOTE]
>The --build command of Kraken2 can take quite a while; it's recommended to specify extra threads with the `--threads` option. Note that on Homebrew installations of Kraken2, using multiple threads appears to be broken.

---
If you want further detail on how to make a custom Kraken2 database, see the manual on the [Kraken2 GitHub](https://github.com/DerrickWood/kraken2)
