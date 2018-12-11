# TPC space-charge distortion conventions



## Distortion of electrons


### Sign conventions
```
distortion = distorted_coordinate - original_coordinate

dr = r_distorted - r_original
drphi = (phi_distorted - phi_original) * r_original
dz = z_distorted - z_original
```
In case of the expected RUN 3 space-charge scenario, electrons are distorted towards the center of the drift volume on each side of the TPC, respectively
* Radial distortions dr
  * Positive sign at inner part of the TPC, negative sign at outer part
* Azimuthal distortions drphi
  * Sign depends on the polarity of the magnetic field (c1 = omega*tau/(1+(omega*tau)^2; omega*tau = -10 * B * v_0 / |E|)
* Distortions along drift direction dz
  * Local
    * A side: positive at small z (< 125 cm), negative at large z (> 125 cm)
    * C side: negative at small abs(z) (> -125 cm), positive at large abs(z) (< -125 cm)
  * Global
    * A side: negative around z = 125 cm, approximates 0 towards the central electrode (z = 0 cm) and readout chamber (z = 250 cm)
    * C side: positive around z = -125 cm, approximates 0 towards the central electrode (z = 0 cm) and readout chamber (z = -250 cm)


### Local distortions
#### Solution of Langevin equation integrated over the length of one z bin towards readout chambers
```
dr = c0 * int(Er/Ez)dz + c1 * int(Ephi/Ez)dz
drphi = c0 * int(Ephi/Ez)dz - c1 * int(Er/Ez)dz
dz = v'(E)/v0 * int(Ez-Ez_0)dz

Ez_0 = -400 V/cm (nominal drift field in z-direction)
Ez = Ez_0 + Ez_spacecharge
v'/v0 = 0.0024 (relative change of the drift velocity)
```
#### Code
The sign of the distortions is determined by the charge of the drifting particle and thus by the direction of the integration, which is the direction towards which the particle drifts. In the numerical integration procedure used, the direction of the integration is given by the sign of gridSizeZ.
* For Er/Ez and Ephi/Ez, i.e. radial and azimuthal distortions, the direction of integration on A and C side is opposite. However, as also the sign of the Ez_0 field is opposite, these two sign changes cancel and the same numerical formula can be used for A and C side.
* In z direction, the direction of integration on A and C side is opposite. As the sign of the relative drift velocity dependence AliTPCPoissonSolver::fgkdvdE depends on the sign of Ez_0, it is also opposite on A and C side and these two sign changes cancel so that the same numerical formula can be used for A and C side.
As the direction of integration (gridSizeZ) is assumed to be positive, which is the case on the A side, the Ez_0 and therefore also AliTPCPoissonSolver::fgkdvdE need to be negative.
##### Local integrals
```
const Double_t ezField = (AliTPCPoissonSolver::fgkCathodeV - AliTPCPoissonSolver::fgkGG) / AliTPCPoissonSolver::fgkTPCZ0;  // nominal drift field in z-direction, -400 V/cm

localIntErOverEz = (gridSizeZ * 0.5) * ((*eR)(i, j) + (*eR)(i, j + 1)) / (ezField + (*eZ)(i, j));
localIntEPhiOverEz = (gridSizeZ * 0.5) * ((*ePhi)(i, j) + (*ePhi)(i, j + 1)) / (ezField + (*eZ)(i, j));
localIntDeltaEz = (gridSizeZ * 0.5) * ((*eZ)(i, j) + (*eZ)(i, j + 1));
```
##### Distortions:
```
const Double_t AliTPCPoissonSolver::fgkdvdE = 0.0024; ///< [cm/V] drift velocity dependency on the E field (from Magboltz for NeCO2N2 at standard environment)

(*distDrDz)(i, j) = fC0 * localIntErOverEz + fC1 * localIntEPhiOverEz;
(*distDPhiRDz)(i, j) = fC0 * localIntEPhiOverEz - fC1 * localIntErOverEz;
(*distDz)(i, j) = localIntDeltaEz * -1 * AliTPCPoissonSolver::fgkdvdE;
```


### Global distortions
Local distortions integrated along the full drift line of the electrons
```
dr = int(dr)driftline
drphi = int(dphi)driftline * r_original
dz = int(dz)driftline
```



## Correction of electrons


### Sign conventions
```
correction = corrected_coordinate - original_coordinate

dr = r_corrected - r_original
drphi = (phi_corrected - phi_original) * r_original
dz = z_corrected - z_original
```
The sign of the corrections should in general be opposite to that of the distortions.


### Local corrections
```
(*corrDrDz)(i, j + 1) = -1 * (*distDrDz)(i, j);
(*corrDPhiRDz)(i, j + 1) = -1 * (*distDPhiRDz)(i, j);
(*corrDz)(i, j + 1) = -1 * (*distDz)(i, j);
```


### Global corrections
#### Regular lookup tables
```
r0 = r(i);
phi0 = phi(m);
z = z(j+1);

r_corrected = r0 + corrDr(i,j+1);
phi_corrected = phi0 + corrDrphi(i,j+1) / r0;
z_corrected = z + corrDz(i,j+1);

corrDr(i,j) = corrDr(i,j+1) + dr(r_corrected, phi_corrected, z_corrected);
corrDrphi(i,j) = ( corrDrphi(i,j+1) / r0 + drphi(r_corrected, phi_corrected, z_corrected) / r_corrected ) * r0;
corrDz(i,j) = corrDz(i,j+1) + dz(r_corrected, phi_corrected, z_corrected);
```
