# Conformal Prediction Bands for Functional Heart-Rate Data (Samsung Health)

This repository contains our final homework project for **Fundamentals of Statistical Learning (Sapienza University of Rome, A.Y. 2025/2026)**.  
The goal of the project is to move from raw smartwatch measurements to a **functional representation** of heart-rate dynamics, and then quantify uncertainty by building **conformal prediction bands** for *future* heart-rate curves. In practice, we take irregular and partially-missing minute-level heart-rate data from **Samsung Health**, turn each day into a smooth curve (separately for day and night), and finally construct prediction bands that aim to cover a new unseen curve with a chosen confidence level.

If you want the full details (assumptions, plots, diagnostics, and results), you can directly check the **report** and the **scripts** included in this repo.

## What we did (high level)

### 1) Data collection & cleaning
- Exported heart-rate data from **Samsung Health** (smartwatch).
- Focused on the most consistent portion of the dataset (**~13 months**, **Jan 13, 2025 → Jan 31, 2026**) because older data had strong gaps due to device changes. 

### 2) Day vs Night functional windows
To capture different physiological patterns and reduce missing segments caused by charging, we created two fixed 10-hour windows:
- **Day:** 09:30–19:30  
- **Night:** 22:00–08:00  

We then applied missingness constraints (e.g., limit large continuous gaps and overall missing share), obtaining **148 day** and **159 night** windows with ~**96%** average coverage. 

### 3) From time series to functional data
- Interpolated small missing segments and applied **smoothing** to obtain functional curves (discrete + continuous representations). 

### 4) Conformal prediction bands (Projection-based Conformal Bands)
We implemented the method in *Lei, Rinaldo, Wasserman* (projection-based conformal bands):
- Split data into **train / calibration (70/30)**  
- Run **FPCA** on train curves to get low-dimensional scores  
- Fit a **Gaussian Mixture Model (GMM)** with **2 components** (selected via BIC)  
- Use calibration conformity scores to build **component-wise bands**  
- Final band = **union/envelope** of component bands (target **α = 0.1**)   

For night curves we also tested a **log-scale fit** to avoid unrealistic negative lower bounds, then transformed the bands back to bpm. 

## Where to look
- **Full report (recommended):** `final_hw_report.pdf` :(https://github.com/1fiocco/Statistical-learning-project/blob/main/final_hw_report.pdf) 
- **Scripts:** see the scripts in the folder.

## Authors
- Riccardo Pugliese (https://github.com/1fiocco)
- Leonardo Suriano (https://github.com/LeoSuriano)
