# BAMmart
A command line utility for parsing bam files and using transcript or gene ids to query for biomart attributes. Currently a work in progress.

## Using BAMmart:
Given a directory containing one or multiple bam files BAMmart can be used currently to fetch information from Ensembl Biomart.

### Getting Started
If you have never used BAMmart before try using `python3 BAMmart.py --help`. There are currently two different commands that can be used. `python3 BAMmart.py helper --search_term [arg]` can help you find a list of filter names that biomart uses to get information on the transcript ids that you are submitting. 
### Submitting Queries
Once you know the query you want to send you may use the second command to perform the search and save the output to a csv file. For example, if I have several dozen bam files and I want to see if any of my ENST ids have RNACentral IDs, you can use the following command: `python3 BAMmart.py query --root_dir bam_dir --attributes rnacentral --output output.csv`. 
### Caveats
Biomart reccommends that queries be sent in batch sizes of no more than 500 at a time. I have found that using a batch size of 500 still results in the request being declined for being too large. The batch size has been hardcoded to 400 IDs at a time for now but will have a parameter to change this in the future. Due batch submission requirement, BAMmart can be slow for data mining large sets of BAM files. Testing this tool has shown that processing 56 BAM files takes BAMmart around 30 minutes worst case.
### Credit
This tool uses pysam for parsing bam files and pybiomart for establishing a connection to biomart to submit requests. 