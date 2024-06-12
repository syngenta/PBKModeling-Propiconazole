# Git Repository for Propiconazole PBK Modeling

## Description

> **Note** More detailed information is provided in the submitted Propiconazole PBK manuscript.

This project developed multispecies (mouse, rat and human) PBK models in the absence of in vivo PK data for the fungicide propiconazole. These models are used in quantifying the interspecies differences between animal models and humans.

A fit-for-purpose read-across approach was integrated with hierarchical clustering - an unsupervised machine learning algorithm, to bridge the gap in in vivo PK data.

The applicability of the machine-learning based read-across approach was validated using penconazole (source) and epoxiconazole (pseudo-unknown target).

Through the machine-learning based read-across approach, difenoconazole was identified as the most appropriate analog for propiconazole.

A mouse PBK model developed and validated for difenoconazole (source) was re-parameterized and applied to propiconazole (target), with the mode of action of CAR/PXR activation incorporated to simulate the in vivo autoinduction of metabolism.
This work represents a substantial advancement toward understanding internal dosimetry and toxicological effects of conazoles, as well as improving risk assessment of propiconazole within the framework of animal alternative safety assessment strategies.

## Getting Started

> **Note** This project involves coding in R. Please refer to the "1. PBK model" in the Appendix of the submitted Propiconazole PBK manuscript for detailed PBK equations.

- The folder "General code" contains the PBK codes that were shared by all the conazoles included in this project.
- The folder of each conazole contains codes to run for three species, mouse, rat and human.
- The folder of clustering contains codes for hierarchical clustering.

## Authors

- Yaoxing Wu ([yaoxing.wu@syngenta.com](mailto:yaoxing.wu@syngenta.com))
- Gabriel Sinclair ([gabriel.sinclair@syngenta.com](mailto:Gabriel.Sinclair@syngenta.com))
- Raghavendhran Avanasi Narasimhan ([raga.avanasi_narasimhan@syngenta.com](mailto:raga.avanasi_narasimhan@syngenta.com))
- Alison Pecquet ([alison.pecquet@syngenta.com](mailto:Alison.Pecquet@syngenta.com))

## Reference

Physiologically based kinetic (PBK) modeling of propiconazole using a machine learning-enhanced read-across approach for interspecies extrapolation.

The article is available online at https://doi.org/10.1016/j.envint.2024.108804 

## Last update
Jun. 12, 2024
