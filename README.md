# Tau-lepton Mass Peak Reconstruction at the Future Circular Collider
This is an implementation of a reconstruction of the tau-lepton mass peak for the Future Circular Collider electron-positron (FCC-ee). This reconstruction employs 3-3-prong events through the decay $\tau^+\tau^-\to\pi^-\pi^+\pi^-\nu_{\tau}$.

This reconstruction has been developed using the Future Circular Collider common analysis framework. Therefore, the `key4hep` stack is required.

## Quick Start

The `key4hep` stack is available and mounted within the LXPLUS Service at CERN. For users with access to LXPLUS, the environment required can be setup by running:

    source setup.sh
This mounts the Future Circular Collider Software (FCCSW) required for the analysis.

In order to access the pre-generated sample used for the analysis, the user should be subscribed to the `fcc-eos-access` e-group within the FCC collaboration.

## Dataframe Production
Once the FCCSW environment is properly mounted, the production of the dataframe containing all the information (*i.e.* the FCC analysis) can be produced by running:

    fccanalysis run Tau3pi_mcSeeded.py
Since the chosen sample contains 100M events, this analysis may take up to an hour to complete. The result is a ROOT file, stored in the `/result` directory containing all the information about the mass reconstruction and the variables (stored as branches within the ROOT file) required for further analyses.

`Tau3pi_mcSeeded.py` corresponds to an analysis for the IDEA detector concept, where the reconstructed events are matched to the MC events; that is, an MC-seeded analysis. A non-MC-seeded analysis would require further adaptations of the steering file to avoid the MC-matching.

## Project Structure

 - **Steering file.** `Tau3pi_mcSeeded.py` is called the steering file. This is a Python script in charge of reading the original sample/s, looping over the events in the sample/s, and writing the results to the final ROOT file. It structures the analysis, but do not perform calculations at the Python level.
 - **Helper file.** `analyzers_Tau3Pi.h` is called the helper file. This is a C++ header file containing all the additional analyzers (*i.e.* functions) required by the steering file. New functions not already defined in the FCCAnalyses framework can be defined in this file. Then, the analyzers are called by the loop in the steering file.
 - **Set-up file.** `setup.sh` is a bash script used to set up the environment required for the analysis. It employs the 2024-03-10 release of the [FCCAnalyses](https://github.com/HEP-FCC/FCCAnalyses).
 - **`/macros` directory.** Once the ROOT file from the steering file has been produced, the C files in this directory can be employed to complete the mass analysis.
 - **`/further_analysis` directory.** Some studies cannot be directly performed using only the original ROOT file produced by the steering file. Some requires new functions, depending on new ideas for the study, or combinations of different branches of that file. This directory provides a centralised C base to extend the study with new ideas.

## Executing the macros (ROOT commands)

The files in the `/macros` and `/further_analyses` directories are C ROOT macros. Their utility is to extend and complete the mass analysis, mainly by producing and plotting results. In order to be executed, these files require the ROOT framework. The ROOT framework can be accessed by typing:

    root

on the command line. However, most of these file require a rootfile to be opened before being executed. To open a rootfile, just type `root` and the file path. For example:

    root result/p8_ee_Ztautau_ecm91.root

`p8_ee_Ztautau_ecm91.root` is the original file produced by the steering file of this project.

Once the rootfile is loaded to ROOT, a macro can be called. All the macros inside `/macros` take an argument, excep for `combHist.C`, which takes no argument. In this project, that argument is always `events`, which is the TTree containing all the branches with the physical information. For example, to call the macro which performs the mass fit:

    root result/p8_ee_Ztautau_ecm91.root
    .x /macros/fit/massFit.C(events)
If a macro is modified or need to be executed again, ROOT should be restarted. To exit ROOT just type:

    .q
And launch ROOT again when the script has been modified.

Some basic ROOT commands to inspect the resulting rootfile are `Print` and `Show`. With `Print`, the data types of the branches can be inspected:

    root result/p8_ee_Ztautau_ecm91.root
    events->Print();
With `Show`, the entries of each branch can be inspected for an specified event. For example:

    root result/p8_ee_Ztautau_ecm91.root
    events->Show(1);
Prints all the data stored in the branches of the event No. `1`, which is passed as an argument.

## What can be done?
Beyond the mass peak reconstruction itself, here is a comprehensive list of what can be achieved with the files provided in this project. All should be within reach with slight or no modification at all.

**Plots**
 - Mass (raw data and fit)
 - Lifetime (raw data and fit)
 - Residuals of $\cos(\beta)$
 - Residuals of $M_{3\pi}$
 - Residuals of 3-pion system 3-momentum
 - Residuals of the secondary vertex positions
 - Correlation matrices of parameters from fit

**Further analyses**
 - Mass corrected by subtracting the energy from radiated photons
 - Mass at MC truth level
 - Mass computed using only MC 3-pion system 3-momentum
 - Mass computed using only MC $\cos(\beta)$
 - Mass computed using only MC secondary vertex
 - **Plots**
	 - Beam spot resolutions ($x$, $y$, $z$)
	 - Vertex detector resolutions ($x$, $y$, $z$)
	 - $\sqrt{s}$

#
This is a work in progress.
