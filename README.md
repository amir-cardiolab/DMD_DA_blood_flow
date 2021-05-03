Codes and data used in the following papers:
	Data-Driven Pulsatile Blood Flow Physics with Dynamic Mode Decomposition
Milad Habibi, Scott Dawson, Amirhossein Arzani
 https://www.mdpi.com/2311-5521/5/3/111
	Integrating multi-fidelity blood flow data with reduced-order data assimilation
	Milad Habibi, Roshan M D'Souza, Scott Dawson, Amirhossein Arzani
	https://arxiv.org/abs/2104.01971
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MATLAB codes (all the input files can be found in the data section):
	DMDc: Dynamic mode decomposition with control.
V: Velocity data, size: M×N
U: Controller input data, size: 1×N
where V and U are the input to the DMDc: DMDc (U,V)
	ROM-KF: Reduced-order modeling Kalman filter.
Gr: Ground truth velocity data, size: M×N
Un: Uncertain velocity data, size: M×N
Ex: Experimental velocity data, size: P×N
H: Observation matrix, size: P×M
Co: The diagonal entries of the diagonal covariance matrix: M×1
where Gr, Un, Ex, H, and Co are the input to the ROM-KF: 
ROM-KF(Gr,Un,Ex,H,Co)
	WOM: Womersley Analytical solution
Fr: Pressure gradient wave form Fourier series’ frequency 
Cn: Pressure gradient wave form Fourier series’ coefficients.
where Fr and Cn are the input to the WOM: WOM (Fr,Cn)
Data:
The datasets are provided in the following Google Drive link and should be placed in the same code directory.
https://drive.google.com/drive/folders/1iBMG4K1_GxfqrHYxM4N9CNnlw4ZWgyDd?usp=sharing 
![image](https://user-images.githubusercontent.com/48529626/116927760-f08a9400-ac10-11eb-8c51-8d1c109e994a.png)



