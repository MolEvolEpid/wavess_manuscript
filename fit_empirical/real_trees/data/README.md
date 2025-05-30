25 fasta files. The first field in the file name (separated by '+') is the patient id as from lanl database. Some patients have samples from more than one genomic region, which were written in different files.

The sequence headers have seven fields (separated by '+') and they are: patient_id+hxb2_start+hxb2_end+accession_number+subtype+days-from-infection+days-from-first-sample

The time stamps for the sequences are relative to time from seroconversion so they are not present in the sequence headers. But you can find the time stamps in the csv file attached.
The first column is the prefix of the fasta file. The timestamps (3rd column) are strings joined with a dash. Be careful if you parse the timestamps automatically since for two cases the first timepoint is a negative value

The fourth column you will find the sequences (accession number) that belong to each timepoint such that sequences in the same timepoint are separated by '+' and different timepoints are separated by '-'.

The references and genes included in each fasta file can be found in the excel file.

Founder file contains the founder sequence for each sample that was used for simulations.
