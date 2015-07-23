###About this project:

- **Project date**: 11/14/2013

- **Objective**: Write a program that encodes realistic bacterial protein-coding genes determined by having it follow specific input and output criteria. See [Project_2_Description.pdf](Project_2_Description.pdf) and [InputOutputCriteriaOutline.txt](InputOutputCriteriaOutline.txt) for a comprehensive overview.

- **Outcome**: The demanding input and output criteria provided to be very challenging for our group to write a program that successfully met all the requirements. The final program was able to encode bacterial protein-coding genes that had a basic bacterial promoter including a Pribnow box, a Shine-Dalgarno sequence, a start and stop codon, and an intrinsic terminator. It was also able to read a protein alignment provided in the input and have the most conserved/variable sites be accurately represented in the output sequences. The program also made the effort to keep the percentage identity between any two protein sequences in the output as low as possible, match the relative codon usage for each amino acid to the codon frequency table provided in the input as closely as possible, and space identical codons, or synonymous codons differing only at the third position, as far apart as possible. The one criterion that our program struggled to manage was having the output contain a functional riboswitch such that at room temperature, the Shine-Dalgarno sequence would only be exposed if the input effector sequence were present. 

- **Contribution**: Outlined the main structural code for the program and overlooked the debugging and code revision process for the group. Identified key functions that would be required for the program, delegated team members to write code for them, and explained the algorithms that would be required to understand the main method and the program’s multiple functions. 


###How to use:

1. Have the required input files in the same directory as [runproject.py](runproject.py). See the *example input files* folder for reference.
2. Run [runproject.py](runproject.py) via command-line: ```python runproject.py```
3. Open *output.fasta* to view the results. (Located in the same directory as [runproject.py](runproject.py))

###Notes:
- This program was created and run using Python 2.7.5 in October/November 2013. Current updates to Python may or may not be compatible with the final version of this code.
