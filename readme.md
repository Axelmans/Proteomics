# Introduction
This repository contains code for a computational biology project course at the University of Antwerp.
The topic was Proteomics, and involved analysing a dataset using MSFragger.
Python code was written to further explore and manipulate the results from this program.

The dataset consists of proteomic profiles from cancer patients belonging to diverse ethnic groups.  
The primary goal was to investigate whether significant differences could be observed between these groups.

This project aims to contribute to unbiased research in proteomics, addressing the still rather prevalent issue that biological studies often disproportionately focus on a single ethnic group.
However, biologically relevant differences between ethnic groups do exist, and being aware of them goes a long way towards better and more personalised care.

# Research Questions  
1. Does an Open Search annotate more spectra than a Closed Search?
2. Which PTM's are present inside each individual and at what frequencies?

# Practical
The Python code uses the __pyteomics__ library. An abstract base class named __Analyser__ is defined in __analyser.py__, which gets derived by more concrete analysers for analysing the __.pepXML__ files.

There is a directory with such files for both searches.

Each derived class reads  __.pepXML__ files, but differs in the way they use their information. 

# Question 1
The relevant file for this question is __stats.py__, here is the general approach of its analysis:
1. Search hits are filtered: only matches with a low enough "expect" score (<= 0.05) count as annotations.
2. Decoy matches and target matches are separated.
3. A warning is printed if the calculated False Discovery Rate (FDR) exceeds its threshold (0.05).
4. All calculated stats are stored in a __.json__ file.

This process is executed for both searches. The results should be as follows:
```
{  
    "Closed": {  
        "Spectra": 5001126,  
        "Expect Threshold": 0.05,  
        "Annotations": 1253250,  
        "Annotation Ratio": 0.251,  
        "Decoy Annotations": 10761,  
        "Target Annotations": 1242489,  
        "False Discovery Rate (FDR)": 0.009,  
        "FDR threshold": 0.005  
  }  
}
```
```
{  
    "Open": {  
        "Spectra": 5601923,  
        "Expect Threshold": 0.05,  
        "Annotations": 1838696,  
        "Annotation Ratio": 0.328,  
        "Decoy Annotations": 9820,  
        "Target Annotations": 1828876,  
        "False Discovery Rate (FDR)": 0.005,  
        "FDR threshold": 0.05  
  }  
}
```
Some key points:
- The Open Search nets a significant larger amount of target annotations.
- Both searches respect their FDR-threshold.
- The amount of decoy hits is close, but the FDR isn't because of the difference in target annotations.

# Question 2
For this question, a local database is used. Filling it takes a considerable amount of time (~3 hours).
Each entry of the database stores a tuple: (individual, ethnic_group, ptm_match, frequency).
The combination of the first 3 items is always unique. In case of a conflict, the frequency is incremented.
The motivation after the database is not having to rerun the match analysis every single time.

The relevant files for this question are __matcher.py__ and __plotter.py__, here is the general approach of its analysis:
1. Again, only highly confident search hits are considered in the analysis.
2. Additionally, decoy matches also get ignored in this analysis.
3. Try to match the mass difference to a combination of PTM's.
4. If the match is close enough, add a database entry.

After the database is filled, we can retrieve the information to generate a plot for each individual, where the x-axis contains all matched PTM's and the y-axis dictates the frequency, here is an example:

![02](https://github.com/user-attachments/assets/42338f23-4096-4be9-a07f-4dd767917049)

Here are some observations:
- The frequencies vary between individuals, likely due to differences in confident search hits and close enough PTM combination matches. The overall shape of the plots seems to match, however.
- Not every PTM is as common (e.g. Methylation is common, Sulfation is not).
- The following  6 PTM's seem to vary the most in frequency: Dimethylation, Dioxidation, Disulfideloss, Nitration, Nitrosylation and Oxidation. Their differences between graphs are most visible.
All other PTM's also show at least some variance, but to a lesser extent.

In general, the PTM matches seem to vary by individual, and not by ethnic group.
