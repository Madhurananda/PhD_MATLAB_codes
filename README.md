# Biologically Inspired Speech Coding & Sound Reconstruction

![MATLAB](https://img.shields.io/badge/MATLAB-research-orange)
![Signal Processing](https://img.shields.io/badge/domain-Speech%20%26%20Auditory%20Coding-blue)
![Research Code](https://img.shields.io/badge/type-PhD%20Research-lightgrey)
![Open Science](https://img.shields.io/badge/open--science-archive-success)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

---

## Overview

This repository contains the **MATLAB research codebase developed during the PhD research of Madhurananda Pahar**, focused on:

* biologically inspired speech coding,
* auditory nerve spike representations,
* event-based signal encoding,
* and sound reconstruction from neural spike trains.

The work investigates how speech signals can be encoded and decoded using **neural spike/event representations**, inspired by biological auditory processing mechanisms.

The repository serves as a **research archive** containing experimental implementations, prototype algorithms, and supporting utilities used across multiple studies on spike-based speech representation and reconstruction.

---

## Research Motivation

Traditional speech processing relies on continuous waveform representations.
This project explores an alternative paradigm:

> Represent speech as **discrete neural spike events**, similar to how the auditory system encodes sound.

Key research questions addressed:

* Can speech be represented using spike/event codes?
* How accurately can audio be reconstructed from spike trains?
* What information is preserved or lost in biologically inspired coding?
* How does spike coding compare with classical signal representations?

---

## Repository Structure

```
PhD_MATLAB_codes/
│
├── MP_Lib/
│   ├── AN Spike Construction/              # Auditory nerve spike generation
│   ├── AN Sound Reconstruction/            # Reconstruction algorithms
│   ├── AN & Onset Sound Reconstruction/    # Onset-based reconstruction studies
│   ├── Signal Comparison & Spectrogram/    # Evaluation & visualisation tools
│   ├── Compression/                        # Spike-based compression experiments
│   ├── Reverberation/                      # Acoustic manipulation utilities
│   ├── Converting Sounds/                  # Audio preprocessing utilities
│   ├── Code Bin/                           # Supporting helper scripts
│   └── Testing/                            # Experimental validation scripts
│
└── README.md
```

The **MP_Lib** directory contains modular MATLAB implementations covering the full experimental workflow:

1. Audio preprocessing
2. Spike/event generation
3. Neural coding experiments
4. Sound reconstruction
5. Signal analysis and evaluation

---

## Main Components

### 1. Auditory Nerve Spike Construction

Implements biologically inspired encoding methods that transform acoustic signals into spike trains representing neural firing activity.

Capabilities include:

* onset detection
* spike timing representation
* amplitude-to-event conversion
* neural-inspired temporal coding

---

### 2. Sound Reconstruction from Spike Codes

Algorithms designed to reconstruct audio signals from spike/event representations.

Explored approaches include:

* auditory nerve reconstruction
* onset-based decoding
* hybrid reconstruction strategies
* perceptual signal comparison

---

### 3. Signal Processing Utilities

Supporting modules for:

* spectrogram generation
* signal comparison metrics
* differentiation and modulation experiments
* reverberation simulation
* silence and synthetic signal generation

---

### 4. Experimental Research Framework

The repository reflects an **iterative experimental PhD workflow**, containing:

* prototype algorithms
* exploratory analyses
* validation experiments
* discussion-oriented implementations

These scripts collectively supported multiple research publications on spike-based speech coding.

---

## Requirements

* MATLAB (recommended R2016b or later)
* Signal Processing Toolbox (recommended)

No external datasets are included.

Users may run experiments using any `.wav` speech signals.

---

## Usage

Because this repository represents a **research archive**, there is no single entry script.

Typical workflow:

1. Navigate to relevant experiment folder inside `MP_Lib/`
2. Load or prepare an audio signal
3. Run spike construction scripts
4. Apply reconstruction algorithms
5. Analyse reconstructed signals using comparison tools

Example workflow:

```
Audio → Spike Encoding → Event Representation → Sound Reconstruction → Evaluation
```

---

## Research Contributions

This codebase contributed to the development of:

* spike/event-based speech representations
* biologically inspired auditory signal processing
* neural coding approaches to speech reconstruction
* alternative paradigms for speech compression and representation

---

## Citation

If you use this repository or related ideas, please cite:

**PhD Thesis**

Pahar, M., 2016.
*A novel sound reconstruction technique based on a spike code (event) representation.*

---

**Conference Publication**

Pahar, M. and Smith, L. S., 2020.
Coding and Decoding Speech using a Biologically Inspired Coding System.
*IEEE Symposium Series on Computational Intelligence (SSCI)*, Canberra, Australia.

```
@INPROCEEDINGS{pahar2020coding,
  author={Pahar, Madhurananda and Smith, Leslie S.},
  booktitle={2020 IEEE Symposium Series on Computational Intelligence (SSCI)}, 
  title={Coding and Decoding Speech using a Biologically Inspired Coding System}, 
  year={2020},
  pages={3025-3032},
  doi={10.1109/SSCI47803.2020.9308328}
}
```

---

## Reproducibility Note

⚠️ This repository is preserved as a **research archive**.

* Scripts reflect experimental PhD development stages.
* Some files are exploratory or prototype implementations.
* Paths and parameters may require adaptation before execution.

---

## Data Availability

No speech datasets are distributed with this repository.

Users must supply their own audio recordings or publicly available speech datasets.

---

## License

Academic research use only unless otherwise specified.

---

## Acknowledgements

This work was developed during doctoral research on biologically inspired speech processing and neural coding systems.

---

⭐ If this repository contributes to your research, please consider starring the project.
