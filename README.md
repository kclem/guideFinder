# guideFinder

```
usage: pams.py [-h] [--exon_file EXON_FILE] [--output_file OUTPUT_FILE]
               [--exon EXON] [--stop_codons STOP_CODONS] [--pams PAMS]
               [--edit_start EDIT_START] [--edit_end EDIT_END]
               [--editing_recognition_seq EDITING_RECOGNITION_SEQ]
               [--editing_resulting_seq EDITING_RESULTING_SEQ]


optional arguments:
  -h, --help            show this help message and exit
  --exon_file EXON_FILE
                        Exon sequence in a file (default: None)
  --output_file OUTPUT_FILE
                        File to which results are written (default: None)
  --exon EXON           Exon sequence (default: None)
  --stop_codons STOP_CODONS
                        Comma-separated sequence of stop codons (default:
                        TGA,TAA,TAG)
  --pams PAMS           Comma-separated sequence of PAMs (default: NGG,NGA)
  --edit_start EDIT_START
                        start index of editing (relative to PAM start)
                        (default: -17)
  --edit_end EDIT_END   end index of editing (relative to PAM start) (this
                        base is not included in the editing window) (default:
                        -12)
  --editing_recognition_seq EDITING_RECOGNITION_SEQ
                        Recognition sequence for editing (default: C)
  --editing_resulting_seq EDITING_RESULTING_SEQ
                        Resulting sequence for editing (default: T)

```
