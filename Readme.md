# An Extended SIR Mathematical Model for Lumpy Skin Disease and Associated Properties

This report presents an extended Susceptible-Exposed-Infectious-Treated-Recovered (SEITR) mathematical
model designed to analyze the transmission dynamics of the Lumpy Skin Disease (LSD) virus in cattle popula-
tions. The model incorporates key aspects of LSD epidemiology and leverages network theory to represent the
complex interactions within the cattle population. The report discusses the model’s formulation, fundamental
properties, implementation, and experimental results. It also outlines potential avenues for future research.

## Introduction
Lumpy Skin Disease (LSD) is a severe disease that affects cattle, causing significant economic losses due to decreased milk production, weight loss, and damage to hides. The disease is caused by the LSD virus, which is primarily transmitted by blood-feeding insects. Understanding the transmission dynamics of the LSD virus is crucial for developing effective control and prevention strategies. To achieve this goal, mathematical modeling has emerged as a powerful tool. In particular, the Susceptible-Exposed-Infectious-Treated-Recovered (SEITR) model, an extension of the classic SIR model, has been widely used to study various diseases. The SEITR model incorporates an exposed stage, representing individuals who have been infected but are not yet infectious, and a treated stage, representing individuals who have been infected and are receiving treatment. In this report, we present an extended SEITR model designed to analyze the transmission dynamics of the LSD virus in cattle populations. Our model incorporates key aspects of LSD epidemiology, including the exposed and treated stages of the disease. Furthermore, our model leverages network theory to represent the complex interactions within the cattle population, which are crucial for understanding the spread of the LSD virus.

## Background
We implemented the SEITR model and this model's application on networks by using the R language. There are five libraries required for this implementation; igraph for network applications, deSolve for solving differential equations, foreach and doParallel for parallelization and Matrix for calculations. To apply SEITR model, we created different networks. We initialized this experiment with a Erdős-Rényi random network, and moved on to a bigger network with more complicated structure. Here the first network is designed to represent a small farmstead where each member of the network is highly connected to the members of the network, and the second network is aimed to replicate a small county with 4 farms that have very high connection within but very low connection in between. In addition to these networks, we account for the newborn members, naturally deceased members and members deceased due to infection in our implementation. To assign status transition, born and deceased nodes, we followed our SEITR model's parameters.
