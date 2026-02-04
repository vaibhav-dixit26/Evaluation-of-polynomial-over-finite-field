# Evaluation of Polynomials over Finite Fields

This repository provides reference implementations of the efficient algorithms proposed in our paper for evaluating multivariate polynomials over finite fields of the form $\mathbb{F}_{q^n}$.

As a proof of concept, we include implementations for both the general case of arbitrary degree $d$ and specialized cases (in particular, $d=2$), along with experimental verification of the theoretical time and memory complexity bounds.

---

## Repository Structure

The repository consists of one helper file and three main SageMath programs.

### Helper File

- **`Poly_helper.sage`**  
  Contains auxiliary functions and common routines used by all main programs.  
  This file must be available in the same directory as the main files and is loaded internally by them.

---

### Main Files

1. **`Poly_evaQuadratic.sage`**  
   Evaluates a quadratic polynomial using our efficient algorithm over the finite field $\mathbb{F}_{q^n}$.

2. **`Poly_evaQuadratic_with_complexity_check.sage`**  
   Evaluates a quadratic polynomial in $n$ variables over $\mathbb{F}_{q^n}$ using our efficient algorithm.  
   In addition, this implementation experimentally verifies the **time complexity** and **memory complexity** against the theoretical bounds derived in the paper.

3. **`Poly_evaGeneral_with_complexity_check.sage`**  
   Evaluates a polynomial of degree $d$ in $n$ variables over $\mathbb{F}_{q^n}$ using our efficient algorithm.  
   This implementation also verifies the **time complexity** and **memory complexity** with respect to the theoretical bounds presented in the paper.

---

## Requirements

- **SageMath** (tested with recent versions)
- A Unix-like environment (Linux or macOS recommended)

Make sure the `sage` executable is available in your environment.

---

## How to Run

Ensure all `.sage` files are in the same directory.  
Run any of the main files from the terminal as follows:

```bash
sage Poly_evaQuadratic.sage
