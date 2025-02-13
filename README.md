# ME459 HW7: Parachute Landing Monte Carlo Simulation  

## Description  
This repository contains **Monte Carlo simulations for parachute landings**, developed for **ME459 Robotics and Unmanned Systems**. The focus is on evaluating the **sensitivity of transition altitude and optimal run-in angles** using statistical methods. The repository includes Python scripts for **simulating parachute descent trajectories**, analyzing landing deviations, and identifying optimal drop conditions.  

## Files Included  

### **Parachute Landing Monte Carlo Simulation**  
- **File:** ME 459 HW 7_1.py  
- **Topics Covered:**  
  - Monte Carlo simulations of parachute descent  
  - Sensitivity analysis of transition altitude  
  - Wind velocity effects on landing location  
  - Statistical evaluation of landing accuracy  

### **Run-in Angle Optimization Simulation**  
- **File:** ME 459 HW 7_2.py  
- **Topics Covered:**  
  - Effect of different run-in angles on landing deviation  
  - No-go zone avoidance analysis  
  - Statistical analysis of mean distance from target  
  - Probability of successful landings outside restricted zones  

### **Homework Problems and Documentation**  
- **File:** ME 459 HW7 2022.pdf  
  - Description of Monte Carlo parachute simulations  
  - Sensitivity analysis of transition altitudes  
  - Run-in angle optimization problem statement  
  - Example output tables and plots  

## Installation  
Ensure Python and the required libraries are installed before running the scripts.  

### **Required Python Packages**  
- numpy  
- matplotlib  
- scipy  
- statistics  
- tqdm  
- tabulate  

To install the necessary packages, run:  

```pip install numpy matplotlib scipy statistics tqdm tabulate```  

## Usage  
1. Open a terminal or Python environment.  
2. Run the desired script using:  

```python ME 459 HW 7_1.py```  
```python ME 459 HW 7_2.py```  

3. View generated plots and statistical tables for analysis.  

## Example Output  

### **Parachute Landing Monte Carlo Simulation**  
- Mean landing distance from target calculated for different transition altitudes  
- Standard deviation and percentage of landings within target radius  

### **Run-in Angle Optimization**  
- Comparison of landing accuracy for different run-in angles  
- Percentage of landings avoiding the no-go zone  

## Contributions  
This repository is designed for educational and research purposes. Feel free to modify and expand upon the existing implementations.  

## License  
This project is open for educational and research use.  

---  

**Author:** Alexander Dowell  

