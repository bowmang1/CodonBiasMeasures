# CodonBiasMeasures
A simple command-line utility to calculate several codon bias measures for a given gene sequence. Currently, RSCU, CAI, and SCUO are implemented.

To compile, simply use:
gcc project.c -lm

Make sure to include the '-lm' to ensure the math library is included. Without this, CAI and SCUO cannot be calculated.

The program is very simple. Upon execution, the user will first be prompted for a .txt file containting the gene sequence that is to be parsed. Two 
samples are included. The user will then be prompted for a .csv file containing the reference organism. Once again, two samples (which correspond to
the provided sequences) are provided. Here is an example: 
![here](https://i.imgur.com/uCGeejj.png).

After this, simply type 'help' for a list of valid commands, and choose your desired Codon Bias Measures to view.

Here are the results of executions of all three measurers for the first test set:

![this](https://i.imgur.com/HGqvnJ3.png)
![that](https://i.imgur.com/jtNCH2V.png)

And here are the results of the executions of all three measures for the second test set:

![other](https://i.imgur.com/Kn7eHLx.png)
![one](https://i.imgur.com/nV3Vng3.png)
